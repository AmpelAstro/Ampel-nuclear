#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : ampel/nuclear/t3/dropboxIO.py
# License           : BSD-3-Clause
# Author            : ?
# Date              : ?
# Last Modified Date: 29.11.2022
# Last Modified By  : simeon.reusch@desy.de

import dropbox
import os
import tempfile
from typing import Any, Optional, Union
from collections.abc import Generator
from functools import lru_cache
from concurrent.futures import ThreadPoolExecutor
from functools import cached_property

import astropy, astropy.time
import backoff
from requests.exceptions import ConnectionError

from ampel.types import UBson, T3Send
from ampel.view.TransientView import TransientView
from ampel.view.T3Store import T3Store
from ampel.struct.UnitResult import UnitResult
from ampel.abstract.AbsPhotoT3Unit import AbsPhotoT3Unit
from ampel.base.AmpelBaseModel import AmpelBaseModel
from ampel.secret.NamedSecret import NamedSecret

from ampel.log.AmpelLogger import AmpelLogger

handle_disconnects = backoff.on_exception(backoff.expo, ConnectionError, max_time=300)


class DropboxUnit(AbsPhotoT3Unit):

    dropbox_token: NamedSecret[str]
    dryRun: bool = False  #: whether to write (and sometimes read) from the dropbox
    dryRunDir: Optional[
        str
    ] = None  # if dryRun is enabled, you must specify a desired temporary directory
    logger: AmpelLogger
    base_location: str = "/mampel"  #:optional diff. directory for testing
    max_connections: int = 20

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        dbx = None
        if self.dryRun:
            if not self.dryRunDir:
                raise ValueError("You have to specify dryRunDir when using dryRun")

        dbx = dropbox.Dropbox(
            self.dropbox_token.get(),
            session=dropbox.create_session(max_connections=self.max_connections + 1),
        )
        try:
            acct = dbx.users_get_current_account()

            if (acct.email == "mitchell.karmen@nyu.edu") or (
                acct.email == "sjoertvanvelzen@gmail.com"
            ):
                self.dbx = dbx
            else:
                raise ValueError(
                    f"Invalid dropbox account {acct.email}. Please check access token."
                )
        except Exception as exc:
            raise ValueError(
                "Invalid access token. Please check login credentials."
            ) from exc

        if self.dryRun:
            if self.dryRunDir:
                if not os.path.exists(self.dryRunDir):
                    os.makedirs(self.dryRunDir)
                self.tmpdir = self.dryRunDir
            else:
                self.tmpdir = tempfile.mkdtemp()
            self.logger.info(
                f"Running in dryRun mode, saving to temporary directories. These are here: {self.tmpdir}"
            )

        self.stats = {"bytes": 0, "files": 0}
        self._pool = ThreadPoolExecutor(max_workers=self.max_connections)
        self._uploads = {}

        if self.base_location == "/mampel":
            self.today = astropy.time.Time.now()
        else:
            self.today = astropy.time.Time(
                self.base_location[13:23]
            )  # shouldn't change from test to test
            if self.dryRun:
                self.create_folder(self.base_location)
            if self.base_location[-15:] not in self.get_files(
                "/mampel/test"
            ):  # hard coding in test dir, 'exists' function on works for files
                self.create_folder(self.base_location)

    def process(
        self, gen: Generator[TransientView, T3Send, None], t3s: Optional[T3Store] = None
    ) -> Union[UBson, UnitResult]:
        """ """
        # DUMMY FUNCTION
        a = 1

    def done(self):
        self.commit()
        self.logger.info(
            f"Wrote {self.stats['bytes']/2**20:.3f} MB in {self.stats['files']} files"
        )

    @handle_disconnects
    def create_folder(self, path):
        if self.dryRun:
            if path.startswith("/"):
                path = path[1:]
            path = os.path.join(self.tmpdir, path)
            os.makedirs(path, exist_ok=True)
            self.logger.info(f"Create {path}")
        else:
            try:
                self.dbx.files_create_folder_v2(path)
            except dropbox.exceptions.ApiError as exc:
                # ignore conflict with existing folder
                if (
                    isinstance(exc.error, dropbox.files.CreateFolderError)
                    and exc.error.is_path()
                    and exc.error.get_path().is_conflict()
                    and exc.error.get_path().get_conflict().is_folder()
                ):
                    ...
                else:
                    raise

    @handle_disconnects
    def put(self, path, payload):
        if self.dryRun:
            if path.startswith("/"):
                path = path[1:]
            path = os.path.join(self.tmpdir, path)
            if not os.path.isdir(os.path.dirname(path)):
                os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "wb") as f:
                f.write(payload)
            self.logger.info(f"Wrote {len(payload)} bytes to {path}")
        else:
            f = self._pool.submit(
                handle_disconnects(self.dbx.files_upload_session_start),
                payload,
                close=True,
            )
            self._uploads[f] = (path, len(payload))
        self.stats["bytes"] += len(payload)
        self.stats["files"] += 1

    def commit(self) -> None:
        entries = []
        for future, (path, offset) in self._uploads.items():
            start_result = future.result()
            entries.append(
                dropbox.files.UploadSessionFinishArg(
                    cursor=dropbox.files.UploadSessionCursor(
                        session_id=start_result.session_id, offset=offset
                    ),
                    commit=dropbox.files.CommitInfo(
                        path=path, mode=dropbox.files.WriteMode.overwrite
                    ),
                )
            )
        if not entries:
            return

        self.logger.info(f"committing {len(entries)} uploads")
        launch_result = handle_disconnects(self.dbx.files_upload_session_finish_batch)(
            entries
        )
        if not launch_result.is_complete():
            status = backoff.on_predicate(
                backoff.expo,
                predicate=lambda job_status: job_status.is_in_progress(),
                max_time=300,
                max_value=5,
            )(handle_disconnects(self.dbx.files_upload_session_finish_batch_check))(
                launch_result.get_async_job_id()
            )
            assert status.is_complete()
            for future, result in zip(
                list(self._uploads.keys()), status.get_complete().entries
            ):
                if result.is_success():
                    del self._uploads[future]
                else:
                    raise RuntimeError(str(result.get_failure()))

    @lru_cache(maxsize=1024)
    @handle_disconnects
    def get_files(self, path):
        if not self.dbx:
            self.dbx = dbx
        entries = []
        if self.dryRun:
            if path.startswith("/"):
                temp_path = path[1:]
            temp_path = os.path.join(self.tmpdir, temp_path)
            try:
                entries = os.listdir(temp_path)
                return entries
            except Exception as e:
                self.logger.info("Directory not found locally, checking dropbox...")

        try:
            list_folder = self.dbx.files_list_folder(path)
        except dropbox.exceptions.ApiError as exc:
            if (
                isinstance(exc.error, dropbox.files.ListFolderError)
                and exc.error.is_path()
                and exc.error.get_path().is_not_folder()
            ):
                raise NotADirectoryError(f"{path} is not a directory") from exc
            else:
                print(path)
                raise

        entries += [entry.name for entry in list_folder.entries]
        has_more = list_folder.has_more
        while has_more:
            list_folder = self.dbx.files_list_folder_continue(list_folder.cursor)
            entries += [entry.name for entry in list_folder.entries]
            has_more = list_folder.has_more

        return entries

    @lru_cache(maxsize=1024)
    @handle_disconnects
    def exists(self, path):
        try:
            self.dbx.files_get_metadata(path)
            return True
        except:
            return False

    @handle_disconnects
    def read_file(self, path):
        if self.dryRun:
            filepath_local = os.path.join(
                self.tmpdir, path[1:] if path.startswith("/") else path
            )

            if os.path.exists(filepath_local):
                with open(filepath_local) as f:
                    return f.read()
            else:
                self.logger.info(
                    "File {} not found locally, searching dropbox".format(
                        filepath_local
                    )
                )

        try:
            return self.dbx.files_download(path)
        except dropbox.exceptions.ApiError as exc:
            if (
                isinstance(exc.error, dropbox.files.DownloadError)
                and exc.error.is_path()
                and exc.error.get_path().is_not_found()
            ):
                raise FileNotFoundError(f"{path} is not a file") from exc
            else:
                raise

#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the functions of utils.py dealing with logging
"""

import PanACoTA.utils as utils
import logging
import os
import pytest
import shutil

GENEPATH = os.path.join("test", "data", "generated-by-utils-tests")

@pytest.fixture(autouse=True)
def setup_teardown_module():
    """
    Remove log files at the end of this test module
    """
    # Init logger to level detail (15)
    # utils.init_logger(LOGFILE_BASE, logging.DEBUG, 'test_utils', verbose=1)
    if os.path.isdir(GENEPATH):
        content = os.listdir(GENEPATH)
        for f in content:
            assert f.startswith(".fuse")
    else:
        os.mkdir(GENEPATH)
    print("setup")

    yield
    shutil.rmtree(GENEPATH, ignore_errors=True)
    print("teardown")


def test_class_filter():
    """
    Check that for a class LessThanFilter(warning), info and debug are allowed,
    but warning, error and critical are not.
    """
    a = utils.LessThanFilter(logging.WARNING)
    record = logging.LogRecord("root", logging.DEBUG, "path", 10, "message", "args", "exc_info")
    assert a.filter(record)
    record = logging.LogRecord("root", logging.INFO, "path", 10, "message", "args", "exc_info")
    assert a.filter(record)
    record = logging.LogRecord("root", logging.WARNING, "path", 10, "message", "args", "exc_info")
    assert not a.filter(record)
    record = logging.LogRecord("root", logging.ERROR, "path", 10, "message", "args", "exc_info")
    assert not a.filter(record)
    record = logging.LogRecord("root", logging.CRITICAL, "path", 10, "message", "args", "exc_info")
    assert not a.filter(record)


def test_class_nolevel_filter():
    """
    Check that for a class LessThanFilter(warning), info and debug are allowed,
    but warning, error and critical are not.
    """
    a = utils.NoLevelFilter(logging.WARNING)
    record = logging.LogRecord("root", logging.DEBUG, "path", 10, "message", "args", "exc_info")
    assert a.filter(record)
    record = logging.LogRecord("root", logging.INFO, "path", 10, "message", "args", "exc_info")
    assert a.filter(record)
    record = logging.LogRecord("root", logging.WARNING, "path", 10, "message", "args", "exc_info")
    assert not a.filter(record)
    record = logging.LogRecord("root", logging.ERROR, "path", 10, "message", "args", "exc_info")
    assert a.filter(record)
    record = logging.LogRecord("root", logging.CRITICAL, "path", 10, "message", "args", "exc_info")
    assert a.filter(record)


def test_logger_default(capsys):
    """
    Test that logger is initialized as expected.
    """
    logfile = os.path.join(GENEPATH, "logfile_test.txt")
    level = logging.DEBUG
    utils.init_logger(logfile, level, "default")
    logger = logging.getLogger("default")
    logger.debug("info debug")
    logger.details("info details")
    logger.info("info info")
    logger.warning("info warning")
    logger.error("info error")
    logger.critical("info critical")
    out, err = capsys.readouterr()
    assert "info info" in out
    assert "info error" in err
    assert "info critical" in err
    with open(logfile + ".log", "r") as logf:
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    with open(logfile + ".log.details") as logf:
        assert logf.readline().endswith(" :: DETAIL :: info details\n")
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    with open(logfile + ".log.debug") as logf:
        assert logf.readline().endswith(" :: DEBUG (from default logger) :: info debug\n")
        assert logf.readline().endswith(" :: DETAIL (from default logger) :: info details\n")
        assert logf.readline().endswith(" :: INFO (from default logger) :: info info\n")
        assert logf.readline().endswith(" :: WARNING (from default logger) :: info warning\n")
        assert logf.readline().endswith(" :: ERROR (from default logger) :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL (from default logger) :: info critical\n")
    with open(logfile + ".log.err", "r") as logf:
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")


def test_logger_verbose1(capsys):
    """
    Test that logger is initialized as expected.
    """
    logfile = os.path.join(GENEPATH, "logfile_test.txt")
    level = logging.DEBUG
    utils.init_logger(logfile, level, "toto", verbose=1)
    logger = logging.getLogger("toto")
    logger.debug("info debug")
    logger.details("info details")
    logger.info("info info")
    logger.warning("info warning")
    logger.error("info error")
    logger.critical("info critical")

    out, err = capsys.readouterr()
    assert "info info" in out
    assert "info warning" in err
    assert "info error" in err
    assert "info critical" in err
    with open(logfile + ".log", "r") as logf:
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    with open(logfile + ".log.details") as logf:
        assert logf.readline().endswith(" :: DETAIL :: info details\n")
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    with open(logfile + ".log.err", "r") as logf:
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    with open(logfile + ".log.debug") as logf:
        assert logf.readline().endswith(" :: DEBUG (from toto logger) :: info debug\n")
        assert logf.readline().endswith(" :: DETAIL (from toto logger) :: info details\n")
        assert logf.readline().endswith(" :: INFO (from toto logger) :: info info\n")
        assert logf.readline().endswith(" :: WARNING (from toto logger) :: info warning\n")
        assert logf.readline().endswith(" :: ERROR (from toto logger) :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL (from toto logger) :: info critical\n")


def test_logger_verbose2(capsys):
    """
    Test that logger is initialized as expected.
    """
    logfile = os.path.join(GENEPATH, "logfile_test.txt")
    level = logging.DEBUG
    utils.init_logger(logfile, level, "toto", verbose=2)
    logger = logging.getLogger("toto")
    logger.debug("info debug")
    logger.details("info details")
    logger.info("info info")
    logger.warning("info warning")
    logger.error("info error")
    logger.critical("info critical")
    out, err = capsys.readouterr()
    assert "info debug" in out
    assert "info details" in out
    assert "info info" in out
    assert "info warning" in err
    assert "info error" in err
    assert "info critical" in err
    with open(logfile + ".log", "r") as logf:
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    with open(logfile + ".log.details") as logf:
        assert logf.readline().endswith(" :: DETAIL :: info details\n")
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    with open(logfile + ".log.err", "r") as logf:
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    with open(logfile + ".log.debug") as logf:
        assert logf.readline().endswith(" :: DEBUG (from toto logger) :: info debug\n")
        assert logf.readline().endswith(" :: DETAIL (from toto logger) :: info details\n")
        assert logf.readline().endswith(" :: INFO (from toto logger) :: info info\n")
        assert logf.readline().endswith(" :: WARNING (from toto logger) :: info warning\n")
        assert logf.readline().endswith(" :: ERROR (from toto logger) :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL (from toto logger) :: info critical\n")


def test_logger_quiet(capsys):
    """
    Test that logger is initialized as expected.
    """
    logfile = os.path.join(GENEPATH, "logfile_test.txt")
    level = logging.DEBUG
    utils.init_logger(logfile, level, "quiet", quiet=True)
    logger = logging.getLogger("quiet")
    logger.debug("info debug")
    logger.details("info details")
    logger.info("info info")
    logger.warning("info warning")
    logger.error("info error")
    logger.critical("info critical")
    out, err = capsys.readouterr()
    print(out)
    print(err)
    assert not out
    assert not err
    with open(logfile + ".log", "r") as logf:
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    with open(logfile + ".log.details") as logf:
        assert logf.readline().endswith(" :: DETAIL :: info details\n")
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    with open(logfile + ".log.err", "r") as logf:
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    with open(logfile + ".log.debug") as logf:
        assert logf.readline().endswith(" :: DEBUG (from quiet logger) :: info debug\n")
        assert logf.readline().endswith(" :: DETAIL (from quiet logger) :: info details\n")
        assert logf.readline().endswith(" :: INFO (from quiet logger) :: info info\n")
        assert logf.readline().endswith(" :: WARNING (from quiet logger) :: info warning\n")
        assert logf.readline().endswith(" :: ERROR (from quiet logger) :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL (from quiet logger) :: info critical\n")


def test_logger_info(capsys):
    """
    Test that when logger is initialized with "INFO" level, it does not return DEBUG info.
    """
    logfile = os.path.join(GENEPATH, "logfile_test.txt")
    level = logging.INFO
    utils.init_logger(logfile, level, "info")
    logger = logging.getLogger("info")
    logger.debug("info debug")
    logger.details("info details")
    logger.info("info info")
    logger.warning("info warning")
    logger.error("info error")
    logger.critical("info critical")
    out, err = capsys.readouterr()
    assert "info info" in out
    assert "info error" in err
    assert "info critical" in err
    with open(logfile + ".log", "r") as logf:
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    assert not os.path.isfile(logfile + ".log.details")
    assert not os.path.isfile(logfile + ".log.debug")
    with open(logfile + ".log.err", "r") as logf:
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")


def test_logger_info_verbose1(capsys):
    """
    Test that when logger is initialized with "INFO" level, it does not return DEBUG info.
    """
    logfile = os.path.join(GENEPATH, "logfile_test.txt")
    level = logging.INFO
    utils.init_logger(logfile, level, "info", verbose=1)
    logger = logging.getLogger("info")
    logger.debug("info debug")
    logger.details("info details")
    logger.info("info info")
    logger.warning("info warning")
    logger.error("info error")
    logger.critical("info critical")
    out, err = capsys.readouterr()
    assert "info info" in out
    assert "info warning" in err
    assert "info error" in err
    assert "info critical" in err
    files = os.listdir(GENEPATH)
    files = [f for f in files if "fuse" not in f]
    assert len(files) == 2
    with open(logfile + ".log", "r") as logf:
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    assert not os.path.isfile(logfile + ".log.details")
    assert not os.path.isfile(logfile + ".log.debug")
    with open(logfile + ".log.err", "r") as logf:
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")


def test_logger_info_verbose2(capsys):
    """
    Test that when logger is initialized with "INFO" level, it does not return DEBUG info.
    """
    logfile = os.path.join(GENEPATH, "logfile_test.txt")
    level = logging.INFO
    utils.init_logger(logfile, level, "info", verbose=2)
    logger = logging.getLogger("info")
    logger.debug("info debug")
    logger.details("info details")
    logger.info("info info")
    logger.warning("info warning")
    logger.error("info error")
    logger.critical("info critical")
    out, err = capsys.readouterr()
    assert "info details" in out
    assert "info info" in out
    assert "info warning" in err
    assert "info error" in err
    assert "info critical" in err
    files = os.listdir(GENEPATH)
    files = [f for f in files if "fuse" not in f]
    assert len(files) == 2
    with open(logfile + ".log", "r") as logf:
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    assert not os.path.isfile(logfile + ".log.details")
    with open(logfile + ".log.err", "r") as logf:
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")


def test_logger_warning(capsys):
    """
    Test that when logger is initialized with "WARNING" level, it does not return
    anything in stdout, as DEBUG and INFO are not returned.
    """
    logfile = os.path.join(GENEPATH, "logfile_test.txt")
    level = logging.WARNING
    utils.init_logger(logfile, level, "warn", verbose=1)
    logger = logging.getLogger("warn")
    logger.debug("info debug")
    logger.details("info details")
    logger.info("info info")
    logger.warning("info warning")
    logger.error("info error")
    logger.critical("info critical")
    out, err = capsys.readouterr()
    assert "info info" in out
    assert "info error" in err
    assert "info critical" in err
    files = os.listdir(GENEPATH)
    files = [f for f in files if "fuse" not in f]
    assert len(files) == 2
    with open(logfile + ".log", "r") as logf:
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    assert not os.path.isfile(logfile + ".log.details")
    with open(logfile + ".log.err", "r") as logf:
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")


def test_logger_warning_verbose1(capsys):
    """
    Test that when logger is initialized with "WARNING" level, it does not return
    anything in stdout, as DEBUG and INFO are not returned.
    """
    logfile = os.path.join(GENEPATH, "logfile_test.txt")
    level = logging.WARNING
    utils.init_logger(logfile, level, "warn", verbose=1)
    logger = logging.getLogger("warn")
    logger.debug("info debug")
    logger.details("info details")
    logger.info("info info")
    logger.warning("info warning")
    logger.error("info error")
    logger.critical("info critical")
    out, err = capsys.readouterr()
    assert "info info" in out
    assert "info error" in err
    assert "info warning" in err
    assert "info critical" in err
    files = os.listdir(GENEPATH)
    files = [f for f in files if "fuse" not in f]
    assert len(files) == 2
    with open(logfile + ".log", "r") as logf:
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    assert not os.path.isfile(logfile + ".log.details")
    with open(logfile + ".log.err", "r") as logf:
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")


def test_logger_warning_verbose2(capsys):
    """
    Test that when logger is initialized with "WARNING" level, it does not return
    anything in stdout, as DEBUG and INFO are not returned.
    """
    logfile = os.path.join(GENEPATH, "logfile_test.txt")
    level = logging.WARNING
    utils.init_logger(logfile, level, "warn", verbose=2)
    logger = logging.getLogger("warn")
    logger.debug("info debug")
    logger.details("info details")
    logger.info("info info")
    logger.warning("info warning")
    logger.error("info error")
    logger.critical("info critical")
    out, err = capsys.readouterr()
    assert "info info" in out
    assert "info details" in out
    assert "info error" in err
    assert "info warning" in err
    assert "info critical" in err
    files = os.listdir(GENEPATH)
    files = [f for f in files if "fuse" not in f]
    assert len(files) == 2
    with open(logfile + ".log", "r") as logf:
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    assert not os.path.isfile(logfile + ".log.details")
    with open(logfile + ".log.err", "r") as logf:
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")


def test_logger_critical(capsys):
    """
    Test that when logger is initialized with "CRITICAL" level, it only returns
    CRITICAL information.

    """
    logfile = os.path.join(GENEPATH, "logfile_test.txt")
    level = logging.CRITICAL
    utils.init_logger(logfile, level, "crit")
    logger = logging.getLogger("crit")
    logger.debug("info debug")
    logger.details("info details")
    logger.info("info info")
    logger.warning("info warning")
    logger.error("info error")
    logger.critical("info critical")
    out, err = capsys.readouterr()
    assert "info info" in out
    assert "info error" in err
    assert "info critical" in err
    files = os.listdir(GENEPATH)
    files = [f for f in files if "fuse" not in f]
    assert len(files) == 2
    with open(logfile + ".log", "r") as logf:
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    assert not os.path.isfile(logfile + ".log.details")
    with open(logfile + ".log.err", "r") as logf:
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")


def test_logger_exists(capsys):
    """
    Test that when the logfiles already exist, it creates new ones with a timestamp added
    """
    logfile = os.path.join(GENEPATH, "logfile_test.txt")
    open(logfile + ".log", "w").close()
    open(logfile + ".log.details", "w").close()
    open(logfile + ".log.debug", "w").close()
    open(logfile + ".log.err", "w").close()
    level = logging.DEBUG
    utils.init_logger(logfile, level, "already_exist", verbose=1)
    logger = logging.getLogger("already_exist")
    logger.debug("info debug")
    logger.details("info details")
    logger.info("info info")
    logger.warning("info warning")
    logger.error("info error")
    logger.critical("info critical")
    out, err = capsys.readouterr()
    assert "info info" in out
    assert "info warning" in err
    assert "info error" in err
    assert "info critical" in err
    # Check that initial log files are still empty
    with open(logfile + ".log", "r") as logf:
        assert logf.readlines() == []
    with open(logfile + ".log.debug", "r") as logf:
        assert logf.readlines() == []
    with open(logfile + ".log.err", "r") as logf:
        assert logf.readlines() == []
    with open(logfile + ".log.details", "r") as logf:
        assert logf.readlines() == []
    # Check for new .log file, remove the one which is empty
    import glob
    logs = glob.glob(logfile + "*" + ".log")
    assert len(logs) == 2
    logs.remove(logfile + ".log")
    with open(logs[0], "r") as logf:
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    # Same thing for .log.err file
    logs_err = glob.glob(logfile + "*" + ".log.err")
    assert len(logs_err) == 2
    logs_err.remove(logfile + ".log.err")
    with open(logs_err[0], "r") as logf:
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")


def test_log_listen(capsys):
    """
    Check that when we log to a queue listener, and then handle the logs
    via logger_thread, the logs appear.
    """
    import multiprocessing
    import threading

    # Create Queue, QueueHandler, and log messages to it
    m = multiprocessing.Manager()
    q = m.Queue()
    qh = logging.handlers.QueueHandler(q)
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    root.handlers = []
    logging.addLevelName(utils.detail_lvl(), "DETAIL")
    root.addHandler(qh)
    logger = logging.getLogger('process')
    logger.debug("debug message")
    logger.log(utils.detail_lvl(), "detail message")
    logger.info("info message")
    logger.warning("warning message")
    logger.error("error message")
    logger.critical("critical message")
    q.put(None)

    # Initialize real logger
    logfile = os.path.join(GENEPATH, "logfile_test.txt")
    utils.init_logger(logfile, 0, '')

    # Listen to QueueHandler and handle messages to stdout/stderr/files
    lp = threading.Thread(target=utils.logger_thread, args=(q,))
    lp.start()
    lp.join()

    out, err = capsys.readouterr()
    assert "info message" in out
    assert "error message" in err
    assert "critical message" in err
    with open(logfile + ".log", "r") as logf:
        assert logf.readline().endswith(" :: INFO :: info message\n")
        assert logf.readline().endswith(" :: WARNING :: warning message\n")
        assert logf.readline().endswith(" :: ERROR :: error message\n")
        assert logf.readline().endswith(" :: CRITICAL :: critical message\n")
    with open(logfile + ".log.details") as logf:
        assert logf.readline().endswith(" :: DETAIL :: detail message\n")
        assert logf.readline().endswith(" :: INFO :: info message\n")
        assert logf.readline().endswith(" :: WARNING :: warning message\n")
        assert logf.readline().endswith(" :: ERROR :: error message\n")
        assert logf.readline().endswith(" :: CRITICAL :: critical message\n")
    with open(logfile + ".log.err", "r") as logf:
        assert logf.readline().endswith(" :: WARNING :: warning message\n")
        assert logf.readline().endswith(" :: ERROR :: error message\n")
        assert logf.readline().endswith(" :: CRITICAL :: critical message\n")


def test_log_no_listen(capsys):
    """
    Check that when we log to a queue listener, but never listen to the queue,
    there is nothing in stderr/stdout/files
    """
    import multiprocessing

    # Create Queue, QueueHandler, and log messages to it
    m = multiprocessing.Manager()
    q = m.Queue()
    qh = logging.handlers.QueueHandler(q)
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    root.handlers = []
    logging.addLevelName(utils.detail_lvl(), "DETAIL")
    root.addHandler(qh)
    logger = logging.getLogger('process')
    logger.debug("debug message")
    logger.log(utils.detail_lvl(), "detail message")
    logger.info("info message")
    logger.warning("warning message")
    logger.error("error message")
    logger.critical("critical message")
    q.put(None)

    # Initialize real logger
    logfile = os.path.join(GENEPATH, "test_log_listen")
    utils.init_logger(logfile, 0, '')

    assert q.qsize() == 7
    out, err = capsys.readouterr()
    assert out == ""
    assert err == ""
    with open(logfile + ".log", "r") as logf:
        assert logf.readlines() == []
    with open(logfile + ".log.details") as logf:
        assert logf.readlines() == []
    with open(logfile + ".log.err", "r") as logf:
        assert logf.readlines() == []


# def test_logger_thread(capsys):
#     """
#     Test that, when we put some fake logs in the Queue given to logger_thread,
#     those logs are given to logger, and printed to err.
#     """
#     import multiprocessing
#     import threading
#     m = multiprocessing.Manager()
#     q = m.Queue()
#     lp = threading.Thread(target=utils.logger_thread, args=(q,))
#     lp.start()
#     q.put(FakeLog("myname", "hello!!"))
#     q.put(FakeLog("other name", "that's me!!!"))
#     q.put(None)
#     lp.join()
#     out, err = capsys.readouterr()
#     print(err)
#     assert out == ""
#     assert "hello!!" in err
#     assert "that's me!!!" in err


# class FakeLog:
#     """
#     Class simulating a logger
#     """

#     def __init__(self, name, text, levelno=100):
#         self.name = name
#         self.text = text
#         self.levelno = levelno
#         self.exc_info = ""
#         self.exc_text = ""
#         self.stack_info = ""

#     def getMessage(self):
#         """
#         returns text of log
#         """
#         return self.text

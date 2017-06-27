#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the functions of utils.py dealing with logging
"""

import genomeAPCAT.utils as utils
import pytest
import logging
import os


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
    logfile = "logfile_test.txt"
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
    assert "info debug\n" in out
    assert "info details\n" not in out
    assert "info info\n" in out
    assert "info warning\n" not in err
    assert "info error\n" in err
    assert "info critical\n" in err
    with open(logfile + ".log", "r") as logf:
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    with open(logfile + ".log.details") as logf:
        assert logf.readline().endswith(" :: DEBUG :: info debug\n")
        assert logf.readline().endswith(" :: DETAIL :: info details\n")
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    with open(logfile + ".log.err", "r") as logf:
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    os.remove(logfile + ".log")
    os.remove(logfile + ".log.err")
    os.remove(logfile + ".log.details")


def test_logger_verbose1(capsys):
    """
    Test that logger is initialized as expected.
    """
    logfile = "logfile_test.txt"
    level = logging.DEBUG
    utils.init_logger(logfile, level, "default", verbose=1)
    logger = logging.getLogger("default")
    logger.debug("info debug")
    logger.details("info details")
    logger.info("info info")
    logger.warning("info warning")
    logger.error("info error")
    logger.critical("info critical")
    out, err = capsys.readouterr()
    assert "info debug\n" in out
    assert "info details\n" not in out
    assert "info info\n" in out
    assert "info warning\n" in err
    assert "info error\n" in err
    assert "info critical\n" in err
    with open(logfile + ".log", "r") as logf:
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    with open(logfile + ".log.details") as logf:
        assert logf.readline().endswith(" :: DEBUG :: info debug\n")
        assert logf.readline().endswith(" :: DETAIL :: info details\n")
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    with open(logfile + ".log.err", "r") as logf:
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    os.remove(logfile + ".log")
    os.remove(logfile + ".log.err")
    os.remove(logfile + ".log.details")


def test_logger_verbose2(capsys):
    """
    Test that logger is initialized as expected.
    """
    logfile = "logfile_test.txt"
    level = logging.DEBUG
    utils.init_logger(logfile, level, "default", verbose=2)
    logger = logging.getLogger("default")
    logger.debug("info debug")
    logger.details("info details")
    logger.info("info info")
    logger.warning("info warning")
    logger.error("info error")
    logger.critical("info critical")
    out, err = capsys.readouterr()
    assert "info debug\n" in out
    assert "info details\n" in out
    assert "info info\n" in out
    assert "info warning\n" in err
    assert "info error\n" in err
    assert "info critical\n" in err
    with open(logfile + ".log", "r") as logf:
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    with open(logfile + ".log.details") as logf:
        assert logf.readline().endswith(" :: DEBUG :: info debug\n")
        assert logf.readline().endswith(" :: DETAIL :: info details\n")
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    with open(logfile + ".log.err", "r") as logf:
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    os.remove(logfile + ".log")
    os.remove(logfile + ".log.err")
    os.remove(logfile + ".log.details")


def test_logger_quiet(capsys):
    """
    Test that logger is initialized as expected.
    """
    logfile = "logfile_test.txt"
    level = logging.DEBUG
    utils.init_logger(logfile, level, "default", quiet=True)
    logger = logging.getLogger("default")
    logger.debug("info debug")
    logger.details("info details")
    logger.info("info info")
    logger.warning("info warning")
    logger.error("info error")
    logger.critical("info critical")
    out, err = capsys.readouterr()
    assert "info debug\n" not in out
    assert "info details\n" not in out
    assert "info info\n" not in out
    assert "info warning\n" not in err
    assert "info error\n" not in err
    assert "info critical\n" not in err
    with open(logfile + ".log", "r") as logf:
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    with open(logfile + ".log.details") as logf:
        assert logf.readline().endswith(" :: DEBUG :: info debug\n")
        assert logf.readline().endswith(" :: DETAIL :: info details\n")
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    with open(logfile + ".log.err", "r") as logf:
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    os.remove(logfile + ".log")
    os.remove(logfile + ".log.err")
    os.remove(logfile + ".log.details")


def test_logger_info(capsys):
    """
    Test that when logger is initialized with "INFO" level, it does not return DEBUG info.
    """
    logfile = "logfile_test.txt"
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
    assert "info debug\n" not in out
    assert "info details\n" not in out
    assert "info info\n" in out
    assert "info warning\n" not in err
    assert "info error\n" in err
    assert "info critical\n" in err
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
    os.remove(logfile + ".log")
    os.remove(logfile + ".log.err")


def test_logger_info_verbose1(capsys):
    """
    Test that when logger is initialized with "INFO" level, it does not return DEBUG info.
    """
    logfile = "logfile_test.txt"
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
    assert "info debug\n" not in out
    assert "info details\n" not in out
    assert "info info\n" in out
    assert "info warning\n" in err
    assert "info error\n" in err
    assert "info critical\n" in err
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
    os.remove(logfile + ".log")
    os.remove(logfile + ".log.err")


def test_logger_info_verbose2(capsys):
    """
    Test that when logger is initialized with "INFO" level, it does not return DEBUG info.
    """
    logfile = "logfile_test.txt"
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
    assert "info debug\n" not in out
    assert "info details\n" not in out
    assert "info info\n" in out
    assert "info warning\n" in err
    assert "info error\n" in err
    assert "info critical\n" in err
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
    os.remove(logfile + ".log")
    os.remove(logfile + ".log.err")


def test_logger_warning(capsys):
    """
    Test that when logger is initialized with "WARNING" level, it does not return
    anything in stdout, as DEBUG and INFO are not returned.
    """
    logfile = "logfile_test.txt"
    level = logging.WARNING
    utils.init_logger(logfile, level, "warn")
    logger = logging.getLogger("warn")
    logger.debug("info debug")
    logger.details("info details")
    logger.info("info info")
    logger.warning("info warning")
    logger.error("info error")
    logger.critical("info critical")
    out, err = capsys.readouterr()
    assert "info debug\n" not in out
    assert "info details\n" not in out
    assert "info info\n" not in out
    assert "info warning\n" not in err
    assert "info error\n" in err
    assert "info critical\n" in err
    with open(logfile + ".log", "r") as logf:
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    assert not os.path.isfile(logfile + ".log.details")
    with open(logfile + ".log.err", "r") as logf:
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    os.remove(logfile + ".log")
    os.remove(logfile + ".log.err")


def test_logger_warning_verbose1(capsys):
    """
    Test that when logger is initialized with "WARNING" level, it does not return
    anything in stdout, as DEBUG and INFO are not returned.
    """
    logfile = "logfile_test.txt"
    level = logging.WARNING
    utils.init_logger(logfile, level, "warn", verbose=1)
    check_warning_verbose(logfile, capsys)


def test_logger_warning_verbose2(capsys):
    """
    Test that when logger is initialized with "WARNING" level, it does not return
    anything in stdout, as DEBUG and INFO are not returned.
    """
    logfile = "logfile_test.txt"
    level = logging.WARNING
    utils.init_logger(logfile, level, "warn", verbose=2)
    check_warning_verbose(logfile, capsys)


def check_warning_verbose(logfile, capsys):
    logger = logging.getLogger("warn")
    logger.debug("info debug")
    logger.details("info details")
    logger.info("info info")
    logger.warning("info warning")
    logger.error("info error")
    logger.critical("info critical")
    out, err = capsys.readouterr()
    assert "info debug\n" not in out
    assert "info details\n" not in out
    assert "info info\n" not in out
    assert "info warning\n" in err
    assert "info error\n" in err
    assert "info critical\n" in err
    with open(logfile + ".log", "r") as logf:
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    assert not os.path.isfile(logfile + ".log.details")
    with open(logfile + ".log.err", "r") as logf:
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    os.remove(logfile + ".log")
    os.remove(logfile + ".log.err")


def test_logger_critical(capsys):
    """
    Test that when logger is initialized with "CRITICAL" level, it only returns
    CRITICAL information.

    """
    logfile = "logfile_test.txt"
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
    assert "info debug\n" not in out
    assert "info details\n" not in out
    assert "info info\n" not in out
    assert "info warning\n" not in err
    assert "info error\n" not in err
    assert "info critical\n" in err
    with open(logfile + ".log", "r") as logf:
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    assert not os.path.isfile(logfile + ".log.details")
    with open(logfile + ".log.err", "r") as logf:
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    os.remove(logfile + ".log")
    os.remove(logfile + ".log.err")


def test_logger_exists(capsys):
    """
    Test that when the logfiles already exist, it creates new ones with a timestamp added
    """
    logfile = "logfile_test.txt"
    open(logfile + ".log", "w").close()
    open(logfile + ".log.details", "w").close()
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
    assert "info debug\n" in out
    assert "info details\n" not in out
    assert "info info\n" in out
    assert "info warning\n" in err
    assert "info error\n" in err
    assert "info critical\n" in err
    with open(logfile + ".log", "r") as logf:
        assert logf.readlines() == []
    with open(logfile + ".log.err", "r") as logf:
        assert logf.readlines() == []
    with open(logfile + ".log.details", "r") as logf:
        assert logf.readlines() == []
    import glob
    logs = glob.glob(logfile + "*" + ".log")
    assert len(logs) == 2
    logs.remove(logfile + ".log")
    with open(logs[0], "r") as logf:
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    logs = glob.glob(logfile + "*" + ".log.err")
    assert len(logs) == 2
    logs.remove(logfile + ".log.err")
    with open(logs[0], "r") as logf:
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    logs = glob.glob(logfile + "*" + ".log.details")
    assert len(logs) == 2
    logs.remove(logfile + ".log.details")
    with open(logs[0], "r") as logf:
        assert logf.readline().endswith(" :: DEBUG :: info debug\n")
        assert logf.readline().endswith(" :: DETAIL :: info details\n")
        assert logf.readline().endswith(" :: INFO :: info info\n")
        assert logf.readline().endswith(" :: WARNING :: info warning\n")
        assert logf.readline().endswith(" :: ERROR :: info error\n")
        assert logf.readline().endswith(" :: CRITICAL :: info critical\n")
    for f in glob.glob(logfile + "*"):
        os.remove(f)


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
    logfile = "test_log_listen"
    utils.init_logger(logfile, 0, '')

    # Listen to QueueHandler and handle messages to stdout/stderr/files
    lp = threading.Thread(target=utils.logger_thread, args=(q,))
    lp.start()
    lp.join()

    out, err = capsys.readouterr()
    assert "debug message\n" in out
    assert "detail message\n" not in out
    assert "info message\n" in out
    assert "warning message\n" not in err
    assert "error message\n" in err
    assert "critical message\n" in err
    with open(logfile + ".log", "r") as logf:
        assert logf.readline().endswith(" :: INFO :: info message\n")
        assert logf.readline().endswith(" :: WARNING :: warning message\n")
        assert logf.readline().endswith(" :: ERROR :: error message\n")
        assert logf.readline().endswith(" :: CRITICAL :: critical message\n")
    with open(logfile + ".log.details") as logf:
        assert logf.readline().endswith(" :: DEBUG :: debug message\n")
        assert logf.readline().endswith(" :: DETAIL :: detail message\n")
        assert logf.readline().endswith(" :: INFO :: info message\n")
        assert logf.readline().endswith(" :: WARNING :: warning message\n")
        assert logf.readline().endswith(" :: ERROR :: error message\n")
        assert logf.readline().endswith(" :: CRITICAL :: critical message\n")
    with open(logfile + ".log.err", "r") as logf:
        assert logf.readline().endswith(" :: WARNING :: warning message\n")
        assert logf.readline().endswith(" :: ERROR :: error message\n")
        assert logf.readline().endswith(" :: CRITICAL :: critical message\n")
    os.remove(logfile + ".log")
    os.remove(logfile + ".log.err")
    os.remove(logfile + ".log.details")


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
    logfile = "test_log_listen"
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
    os.remove(logfile + ".log")
    os.remove(logfile + ".log.err")
    os.remove(logfile + ".log.details")


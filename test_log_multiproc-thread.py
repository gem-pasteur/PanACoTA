import logging
import logging.config
import logging.handlers
from multiprocessing import Process, Queue
import random
import threading
import time

import genomeAPCAT.utils as utils

def logger_thread(q):
    while True:
        record = q.get()
        if record is None:
            break
        logger = logging.getLogger(record.name)
        logger.handle(record)


def worker_process(q):
    qh = logging.handlers.QueueHandler(q)
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    root.addHandler(qh)
    levels = [logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR,
              logging.CRITICAL]
    loggers = ['foo', 'foo.bar', 'foo.bar.baz',
               'spam', 'spam.ham', 'spam.ham.eggs']
    for i in range(10):
        lvl = random.choice(levels)
        logger = logging.getLogger('process')
        logger.log(lvl, 'Message no. %d', i)


if __name__ == '__main__':
    q = Queue()
    # d = {
    #     'version': 1,
    #     'formatters': {
    #         'detailed': {
    #             'class': 'logging.Formatter',
    #             'format': '%(asctime)s %(name)-15s %(levelname)-8s %(processName)-10s %(message)s'
    #         }
    #     },
    #     'handlers': {
    #         'console': {
    #             'class': 'logging.StreamHandler',
    #             'level': 'INFO',
    #         },
    #         'file': {
    #             'class': 'logging.FileHandler',
    #             'filename': 'mplog.log',
    #             'mode': 'w',
    #             'formatter': 'detailed',
    #         },
    #         'foofile': {
    #             'class': 'logging.FileHandler',
    #             'filename': 'mplog-foo.log',
    #             'mode': 'w',
    #             'formatter': 'detailed',
    #         },
    #         'errors': {
    #             'class': 'logging.FileHandler',
    #             'filename': 'mplog-errors.log',
    #             'mode': 'w',
    #             'level': 'ERROR',
    #             'formatter': 'detailed',
    #         },
    #     },
    #     'loggers': {
    #         'foo': {
    #             'handlers': ['foofile']
    #         }
    #     },
    #     'root': {
    #         'level': 'DEBUG',
    #         'handlers': ['console', 'file', 'errors']
    #     },
    # }


    workers = []
    for i in range(5):
        wp = Process(target=worker_process, name='worker %d' % (i + 1), args=(q,))
        workers.append(wp)
        wp.start()


    utils.init_logger('mplog', logging.DEBUG, '')
    logger = logging.getLogger('main')
    # logging.config.dictConfig(d)
    lp = threading.Thread(target=logger_thread, args=(q,))
    lp.start()
    logger.info("workers started")
    # At this point, the main process could do some useful work of its own
    # Once it's done that, it can wait for the workers to terminate...
    for wp in workers:
        wp.join()
    # And now tell the logging thread to finish up, too
    logger.info("workers ended")
    q.put(None)
    lp.join()
    logger.details("THE END!")
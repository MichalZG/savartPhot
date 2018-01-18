import logging


def get_logger(name, config_data):
    logger = logging.getLogger(name)
    handler = logging.FileHandler('var/logs/{}.log'.format(name))
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)

    if config_data:
        for key in config_data.keys:
            logger.debug('CONFIG {0}:\t{1}'.format(key, config_data.get(key)))
    return logger

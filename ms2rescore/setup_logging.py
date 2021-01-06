import logging


def setup_logging(passed_level):
    log_mapping = {
        'critical': logging.CRITICAL,
        'error': logging.ERROR,
        'warning': logging.WARNING,
        'info': logging.INFO,
        'debug': logging.DEBUG,
    }

    if passed_level not in log_mapping:
        print(
            "Invalid log level. Should be one of the following: ",
            ', '.join(log_mapping.keys())
        )
        exit(1)

    logging.basicConfig(
        format='%(asctime)s // %(levelname)s // %(name)s // %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        level=log_mapping[passed_level]
    )

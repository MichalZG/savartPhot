from src.utils.logger import get_logger


class PipelineBase(object):
    def __init__(self, log_file_name, config):
        self.logger = self._get_logger(log_file_name, config)

    @staticmethod
    def _get_logger(log_file_name, config):
        if log_file_name:
            return get_logger(log_file_name, config)
        return None

    def log(self, message):
        if self.logger:
            self.logger.log(message)

    def debug(self, message):
        if self.logger:
            self.logger.debug(message)

    def info(self, message):
        if self.logger:
            self.logger.info(message)

    def warning(self, message):
        if self.logger:
            self.logger.warning(message)

    def error(self, message):
        if self.logger:
            self.logger.error(message)

    def critical(self, message):
        if self.logger:
            self.logger.critical(message)

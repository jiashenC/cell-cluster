import time


class Log:
    def __init__(self):
        self.pre = time.time()

    def time(self, status):
        print '%s: %f, %f' % (status, time.time(), (time.time() - self.pre))
        self.pre = time.time()

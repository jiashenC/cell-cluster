import time


class Log:
    def __init__(self):
        self.pre = time.time()

    def time(self, status):
        print '%s: %f, %f' % (status, time.time(), (time.time() - self.pre))
        self.pre = time.time()

    def update(self):
        self.pre = time.time()


def output(step, data):
    print "*****************************************************"
    print
    print step
    print data
    print "*****************************************************"
    print

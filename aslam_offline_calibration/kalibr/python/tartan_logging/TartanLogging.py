class TartanLogger(object):
    def __init__(self,ObservationDatabase,camera,stats):
        self.ObservationDatabase_ = ObservationDatabase
        self.camera_ = camera
        self.stats_ = stats
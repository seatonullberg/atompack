class PositionError(Exception):
    """Base class for position errors."""

    def __init__(self, position):
        self.position = position


class PositionOccupiedError(PositionError):
    """Error raised when a position is unexpectedly occupied."""

    def __init__(self, position):
        super().__init__(position)

    def __str__(self):
        return "Position '{}' is occupied.".format(self.position)


class PositionUnoccupiedError(PositionError):
    """Error raised when a position is unexpectedly unoccupied."""

    def __init__(self, position):
        super().__init__(position)

    def __str__(self):
        return "Position '{}' is unoccupied.".format(self.position)


class PositionOutsideError(PositionError):
    """Error raised when a position is unexpectedly outside of the bounding box."""

    def __init__(self, position, basis):
        self.basis = basis
        super().__init__(position)

    def __str__(self):
        return "Position '{}' is outside of the allowable bounding area '{}'.".format(self.position, self.basis)

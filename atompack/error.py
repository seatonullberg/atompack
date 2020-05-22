class BasePositionError(Exception):
    """The base class for errors which deal with positions."""

    def __init__(self, position, tolerance):
        self.position = position
        self.tolerance = tolerance


class PositionOccupiedError(BasePositionError):
    """The error raised when a position is occupied but should not be."""

    def __init__(self, position, tolerance):
        super().__init__(position, tolerance)

    def __str__(self):
        return "There is an `Atom` within the tolerance radius '{}' of position '{}'.".format(
            self.tolerance, self.position)


class PositionUnoccupiedError(BasePositionError):
    """The error raised when a position is not occupied but should be."""

    def __init__(self, position, tolerance):
        super().__init__(position, tolerance)

    def __str__(self):
        return "There are no `Atom`s within the tolerance radius '{}' of position '{}'.".format(
            self.tolerance, self.position)

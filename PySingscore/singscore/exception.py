class InvalidNormalisation(Exception):
    """
    Raised when an invalid normalisation string is entered
    """
    pass

class InvalidIDType(Exception):
    """
    Raised when the gene identifier types don't match
    """
    pass

class InvalidGrid(Exception):
    """
    Raised when the grid for graphs is not of sufficient dimensions
    """
    pass
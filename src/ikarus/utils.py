import pandas as pd


def intersection_fun(x):
    return pd.concat(x, axis=1, join="inner")


def union_fun(x):
    return pd.concat(x, axis=1, join="outer")

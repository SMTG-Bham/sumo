import math


def multl_const(const, list):
    return [const * item for item in list]


def divl_const(const, list):
    return [item / const for item in list]


def minl_const(const, list):
    return [item - const for item in list]


def multl(list_x, list_y):
    return [x_i * y_i for x_i, y_i in zip(list_x, list_y)]


def divl(list_x, list_y):
    return [x_i / y_i for x_i, y_i in zip(list_x, list_y)]


def minl(list_x, list_y):
    return [x_i - y_i for x_i, y_i in zip(list_x, list_y)]


def addl(list_x, list_y):
    if not list_x:
        return list_y  # useful for adding to an empty list
    return [x_i + y_i for x_i, y_i in zip(list_x, list_y)]


def mag(v):
    return math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)


def ceil_int(x):
    return int(math.ceil(x))


def floor_int(x):
    return int(math.floor(x))


def dot(list_x, list_y):
    return sum([i*j for (i, j) in zip(list_x, list_y)])


def angle(u, v):
    top = dot(u, v)
    bottom = mag(u) * mag(v)
    angle_in_rad = math.acos(top/bottom)
    return math.degrees(angle_in_rad)

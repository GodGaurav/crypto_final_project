"""
Program that demonstrates basic elliptic curve computations
Written by: Qi Ying Lim
"""
INFINITY = (0.5, 0.5)  # represents the additive identity on the elliptic curve
p = 0
b = 0
c = 0
g = (0, 0)


def set_param(b_val, c_val, prime):
    """
    Sets global parameters, where the elliptic curve y^2 = x^3 +bx + c (mod p)
    is defined on the f ield Z_p
    :param b_val: coefficient of x
    :param c_val: coefficient of 1
    :param prime: prime p
    """
    if (4 * (b_val ** 2) + 27 * c_val) % prime == 0:
        raise ValueError("parameters are not valid")
    else:
        global b
        global p
        global c
        b = b_val
        p = prime
        c = c_val


def set_p256_param():
    """
    void method sets parameters of the elliptic curve to P256-standard
    """
    p256_b = -3
    p256_c = 41058363725152142129326129780047268409114441015993725554835256314039467401291
    p256_prime = 115792089210356248762697446949407573530086143415290314195533631308867097853951
    x = 48439561293906451759052585252797914202762949526041747995844080717082404635286
    y = 36134250956749795798585127919587881956611106672985015071877198253568414405109
    global g
    g = (x, y)

    set_param(p256_b, p256_c, p256_prime)


def add(p1, p2):
    """
    Adds two points on an elliptic curve E(Z_p): y^2 = x^3 + bx + c
    :param p1: a tuple representing point 1
    :param p2: a tuple representing point 2
    :return: the sum p1 + p2
    """
    (x_1, y_1) = p1
    (x_2, y_2) = p2

    if additive_inverse(p1, p2):
        p3 = INFINITY
    elif p1 == INFINITY:
        p3 = p2
    elif p2 == INFINITY:
        p3 = p1
    else:
        m = find_m(p1, p2)
        x_3 = (m ** 2 - x_1 - x_2) % p
        y_3 = (m * (x_1 - x_3) - y_1) % p
        p3 = (x_3, y_3)
    return p3

def subtract(p1, p2):
    """
    Subtracts two points on the elliptc curve
    :param p1: tuple representing first point
    :param p2: tuple representing the second point
    :return: p1 - p2
    """
    (x,y) = p2
    neg_p2 = (x, ((-1)*y)%p)
    return add(p1, neg_p2)

def times(n, point):
    """
    Returns product of n and point with Montgomery ladder algorithm
    :param n: integer, where output is n times of point
    :param point: tuple representing point on elliptic curve
    :return: n*point
    """

    # for i in xrange(n - 1):
    #     point = add(point, point)

    if n == 0 or point[0] == 0.5:
        return INFINITY
    elif n == 1:
        return point
    else:
        n_b = bin(n)[2:]
        R_0 = point
        R_1 = add(point, point)
        for char in n_b[1:]:
            if int(char, 2) == 0:
                R_1 = add(R_0, R_1)
                R_0 = add(R_0, R_0)
            else:
                R_0 = add(R_0, R_1)
                R_1 = add(R_1, R_1)
    return R_0


def is_in_graph(g):
    """
    Verifies that point g is on graph of elliptic curve E(Z_p) defined by
    y^2 = x^3 +bx + c (mod p)
    :param g: tuple representing point
    :return: boolean True if point is on graph, False otherwise
    """
    (x, y) = g
    y2 = (y ** 2) % p
    curve = (pow(x, 3, p) + (b * x) + c) % p
    return y2 == curve


def encode(m, K):
    """
    Encodes message m onto the elliptic curve
    :param m: integer representing message
    :param K: integer that determines the failure rate of this encoding
    :return: the encoded form of m
    """

    if (m + 1 * K) >= p:
        raise ValueError('K is too large for given parameters')
    else:
        for j in range(K):
            x = (m * K + j) % p
            n = pow(x, 3, p) + (b * x + c) % p  # computes y^2

            if is_quad_residue(n):  # if n is a quadratic residue modulo p
                y = square_root(n)[0]
                return (x, y)

            else:
                if j == K - 1:  # encoding has failed
                    return "Failed to encode message"


def decode(c, K):
    """
    from point on curve, obtain original plaintext
    :param c: tuple representing a point on the curve
    :param K: integer determining the failure rate of encoding
    :return: original plaintext m
    """
    return c[0]//K


def extended_euclid(num1, num2):
    # taken from Professor Straubing's euclid_etc.py file
    # given values x,y return
    # (gcd(x,y),r,s) where r*x+s*y=gcd(x,y)
    (a, b, aa, bb) = (0, 1, 1, 0)
    while num2 != 0:
        (q, r) = divmod(num1, num2)
        (a, b, aa, bb) = (aa - q * a, bb - q * b, a, b)
        (num1, num2) = (num2, r)
    return (num1, (aa, bb))


def additive_inverse(p1, p2):
    """
    Determines if p1 = -p2,
    :param p1: tuple representing point 1
    :param p2: tuple representing point 2
    :return: True if p1 = -p2, False otherwise
    """
    (x_1, y_1) = p1
    (x_2, y_2) = p2

    if x_1 == x_2 and y_1 == ((-1) * y_2) % p:
        return True
    else:
        return False


def find_m(p1, p2):
    """
    Finds the gradient m
    """
    (x_1, y_1) = p1
    (x_2, y_2) = p2

    if p1 == p2:
        numerator = 3 * (x_1 ** 2) + b
        denominator = extended_euclid(2 * y_1, p)[1][0]
    else:
        numerator = y_2 - y_1
        denominator = extended_euclid(x_2 - x_1, p)[1][0]

    return (numerator * denominator) % p


def is_quad_residue(a):
    """
    Determines if a is a quadratic residue modulo p
    :return: True if a is a quadratic residue
    """
    if a % p == 0:
        return True
    return pow(a, (p - 1) / 2, p) == 1


def square_root(n):
    """
    Uses Shanks Tonelli algorithm to find square root n modulo p
    :param n: integer that is a quadratic residue modulo p
    :return: square-root of n modulo p
    """
    if not is_quad_residue(n):
        return "no square roots exist"
    # express p-1 = (2^e)*s
    s = p - 1
    e = 0

    while s % 2 == 0:
        e += 1
        s /= 2
    q = 1
    while is_quad_residue(q):
        q += 1

    y = pow(n, (s + 1) / 2, p)
    B = pow(n, s, p)
    G = pow(q, s, p)
    r = e

    while B != 1:
        m = 0
        exp = 1
        for m in range(r):
            if pow(b, exp, p) == 1:
                break
            exp *= 2

        if m == 0:
            break
        else:
            y = (y * pow(G, pow(2, r - m - 1), p)) % p;
            G = pow(G, pow(2, r - m), p);
            B = (B * G) % p;
            r = m
    y = min(y, p - y)
    return y, p - y




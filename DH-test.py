from Crypto.PublicKey import ECC
import ecc_simple
"""
Compares output with pycryptodome's output,
Diffie-Hellman demonstration
By: Qi Ying Lim
"""


def format(eccPoint):
    """
    Formats EccPoint object so as to allow comparison with my own ecc_simple function outputs
    :param eccPoint: EccPoint Object
    :return: a tuple of integers
    """
    return (int(eccPoint.x),int(eccPoint.y))

"""
First use pycryptodome's ECC module to generate private keys for Alice and Bob
and compute the shared key with pycryptodome's methods
"""

keyA = ECC.generate(curve='P-256')  # Alices' Key
keyB = ECC.generate(curve='P-256')  # Bob's Key

privateA = keyA.d  # Alice's Private Key
privateB = keyB.d  # Bob's Private Key

publicA = keyA.pointQ  # Alice's Public Key
publicB = keyB.pointQ  # Bob's Public Key

shared_keyA = format(publicB * privateA)  # shared key as computed by each
shared_keyB = format(publicA * privateB)

"""
Verifying ecc_simple methods
"""
ecc_simple.set_p256_param()  # set parameters of the p256 curve
g = ecc_simple.g  # base point

# I am using the same private keys generated by pycryptodome, namely privateA and privateB
publicA_1 = ecc_simple.times(int(privateA), g)  # Alice's Public Key
publicB_1 = ecc_simple.times(int(privateB), g)  # Bob's Public Key

shared_keyA_1 = ecc_simple.times(int(privateA), publicB_1)  # shared key as computed by each
shared_keyB_1 = ecc_simple.times(int(privateA), publicB_1)

print publicA_1
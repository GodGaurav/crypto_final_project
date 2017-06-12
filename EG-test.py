from ecc_simple import *
import random
"""
ElGamal demonstration
By: Qi Ying Lim
"""

set_p256_param()
p = 115792089210356248762697446949407573530086143415290314195533631308867097853951
x = 48439561293906451759052585252797914202762949526041747995844080717082404635286
y = 36134250956749795798585127919587881956611106672985015071877198253568414405109
g = (x, y)

a = random.randrange(1, p)  # Alice's Private Key
b = random.randrange(1, p)   # Bob's Private Key

public_Bob = times(b, g) # Bob's public Key

# Alice wants to send message m to Bob
m = 32454324564335645345643
public_Alice = times(a, g)   # Alice's public Key
shared_key_A = times(a, public_Bob) # shared Key as computed by Alice

c_1 = public_Alice
c_2 = add(encode(m, 100), shared_key_A)
ciphertext = (c_1, c_2)    # Alice sends this to Bob

#To decrypt, Bob computes
shared_key_B = times(b, ciphertext[0])    # shared Key as computed by Bob
plaintext = subtract(ciphertext[1], shared_key_B) # obtains plaintext,
m = decode(plaintext, 100)    # decodes point on curve to retrieve m

print m
import ecc
import curves
import sys
import random

challenge_points = {
    'ecp131':
        [
            int('3AA6F004FC62E2DA1ED0BFB62C3FFB568', 16),
            int('09C21C284BA8A445BB2701BF55E3A67ED', 16)
        ],
    'ecp79':
        [
            int('0679834CEFB7215DC365', 16),
            int('4084BC50388C4E6FDFAB', 16)
        ],
        'ecp89':
        [
            int('00DE1AA94FF94DB64E763E2D', 16),
            int('002A44C4C2D4EE27FA0A4BA9', 16)
        ]
}

def calc_sqrt_val(p):
    if p % 4 == 3:
        print('//(p + 1) // 4')
        return (p + 1) // 4
    elif p % 8 == 5:
        print('//(p - 5) // 8')
        return (p - 5) // 8
    else:
        raise "sqrt error"

def inv(x, p):

    return pow(x, p - 2, p)

def calc_k(p):

    r = (2**160)

    rinv = inv(r % p, p)

    k = (r * rinv - 1) // p

    return k

def get_params(name):

    if name == "ecp131":
        return curves.ecp131
    if name == "ecp79":
        return curves.ecp79
    if name == "ecp89":
        return curves.ecp89
    else:
        raise "Invalid curve"

# Shorter name used some places in code
def get_alt_name(name):

    if name == "ecp131":
        return "p131"
    if name == "ecp79":
        return "p79"
    if name == "ecp89":
        return "p89"
    else:
        raise "Invalid curve"
    
def to_montgomery(n, p):
    r = 2**160

    return (n * r) % p

def to_hex(n):

    a0 = n % (2**64)
    a1 = (n >> 64) % (2**64)
    a2 = (n >> 128)

    return f"{{{{0x{a0:x}, 0x{a1:x}, 0x{a2:x}}}}}"


def randint(n):
    return random.randint(0, 2**256) % n

def main():

    curve_name = sys.argv[1]
    curve = curves.getCurveByName(curve_name)
    params = get_params(curve_name)

    qx = challenge_points[curve_name][0]
    qy = challenge_points[curve_name][1]

    if not curve.verifyPoint(ecc.ECPoint(qx, qy)):
        raise "Error"

    print(f"CurveParameters _{curve_name} = {{")
    print(f"  .p = {to_hex(params.p)},")
    print(f"  .a = {to_hex(to_montgomery(params.a, params.p))},")
    print(f"  .b = {to_hex(to_montgomery(params.b, params.p))},")
    print(f"  .n = {to_hex(params.n)},")
    
    print(f"  .gx = {to_hex(to_montgomery(params.x, params.p))},")
    print(f"  .gy = {to_hex(to_montgomery(params.y, params.p))},")
    
    print(f"  .qx = {to_hex(to_montgomery(qx, params.p))},")
    print(f"  .qy = {to_hex(to_montgomery(qy, params.p))},")
    
    print(f"  .k = {to_hex(calc_k(params.p))},")
    print(f"  .one = {to_hex(to_montgomery(1, params.p))},")
    print(f"  .two = {to_hex(to_montgomery(2, params.p))},")
    
    print(f"  .p_minus_2 = {to_hex(params.p - 2)},")
    print(f"  .sqrt = {to_hex(calc_sqrt_val(params.p))},")
    
    print(f"  .r = {to_hex(to_montgomery(1, params.p))},")
    
    print(f"  .r2 = {to_hex(to_montgomery(2**160, params.p))},")
    
    print(f"  .bits = {params.p.bit_length()},")
    
    print(f"  .words = {(params.p.bit_length() + 63) // 64},")

    print(f"  .name = \"{curve_name}\",")

    print("};")

    alt_name = get_alt_name(curve_name)

    print()
    print(f"__constant__ uint131_t _{alt_name}_p = {to_hex(params.p)};")
    print(f"__constant__ uint131_t _{alt_name}_k = {to_hex(calc_k(params.p))};")
    print(f"__constant__ uint131_t _{alt_name}_r2 = {to_hex(to_montgomery(2**160, params.p))};")
    print(f"__constant__ uint131_t _{alt_name}_one = {to_hex(to_montgomery(1, params.p))};")
    print(f"__constant__ uint131_t _{alt_name}_a = {to_hex(to_montgomery(params.a, params.p))};")
    print(f"__constant__ uint131_t _{alt_name}_b = {to_hex(to_montgomery(params.b, params.p))};")
    print()

    exponents = []
    r_points = []

    random.seed(1)

    p = curve.bp
    q = ecc.ECPoint(qx, qy)
    for _ in range(32):
        a = randint(params.n)
        b = randint(params.n)
        exponents.append((a, b))

        rp = curve.add(curve.multiply(a, p), curve.multiply(b, q))
        r_points.append(rp)

    print()
    print(f"_a_str = [")
    for e in exponents:
        print(f"\"{e[0]:x}\",")
    print("]")
    print()
    print(f"_b_str = [")
    for e in exponents:
        print(f"\"{e[1]:x}\",")
    print("]")
    print()
    print()

    print(f"std::string _{alt_name}_a_str[] = {{")
    for e in exponents:
        print(f"\"{e[0]:x}\",")
    print("};")
    print()
    print(f"std::string _{alt_name}_b_str[] = {{")
    for e in exponents:
        print(f"\"{e[1]:x}\",")
    print("};")
    print()
    print(f"std::string _{alt_name}_x_str[] = {{")
    for r in r_points:
        print(f"\"{to_montgomery(r.x, curve.p):x}\",")
    print("};")
    print() 
    print(f"std::string _{alt_name}_y_str[] = {{")
    for r in r_points:
        print(f"\"{to_montgomery(r.y, curve.p):x}\",")
    print("};")

if __name__ == "__main__":
    main()





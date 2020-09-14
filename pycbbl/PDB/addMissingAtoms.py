import math, warnings
from Bio.PDB.Atom import Atom
from Bio.PDB.vectors import Vector, rotaxis, calc_dihedral, calc_angle


def add_terminal_OXT(chain, C_OXT_length=1.23):
    """Adds a terminal oxygen atom ('OXT') to the last residue of the chain."""

    rad = 180.0 / math.pi

    residues = [r for r in chain.get_residues()]
    last_residue = residues[-1]

    # Check wether the last residues has an OXT atom
    for a in last_residue.get_atoms():
        if a.name == 'OXT':
            raise ValueError('Last residue of input chain already contains an OXT atom')
        if a.name == 'N':
            N_resRef = a
        elif a.name == 'CA':
            CA_resRef = a
        elif a.name == 'C':
            C_resRef = a
        elif a.name == 'O':
            O_resRef = a

    n_vec = N_resRef.get_vector()
    ca_vec = CA_resRef.get_vector()
    c_vec = C_resRef.get_vector()
    o_vec = O_resRef.get_vector()

    # geometry to bring together residue
    CA_C_OXT_angle = calc_angle(ca_vec, c_vec, o_vec) * rad
    N_CA_C_O_diangle = calc_dihedral(n_vec, ca_vec, c_vec, o_vec) * rad
    N_CA_C_OXT_diangle = N_CA_C_O_diangle - 180.0
    if N_CA_C_O_diangle < 0:
        N_CA_C_OXT_diangle = N_CA_C_O_diangle + 180.0

    # OXT atom creation
    OXT_coord = calculateCoordinates(
        N_resRef, CA_resRef, C_resRef, C_OXT_length, CA_C_OXT_angle, N_CA_C_OXT_diangle
    )
    OXT = Atom("OXT", OXT_coord, 0.0, 1.0, " ", "OXT", 0, "O")

    # modify last residue of the structure to contain the OXT atom
    last_residue.add(OXT)

def calculateCoordinates(refA, refB, refC, L, ang, di):
    AV = refA.get_vector()
    BV = refB.get_vector()
    CV = refC.get_vector()

    CA = AV - CV
    CB = BV - CV

    ##CA vector
    AX = CA[0]
    AY = CA[1]
    AZ = CA[2]

    ##CB vector
    BX = CB[0]
    BY = CB[1]
    BZ = CB[2]

    ##Plane Parameters
    A = (AY * BZ) - (AZ * BY)
    B = (AZ * BX) - (AX * BZ)
    G = (AX * BY) - (AY * BX)

    ##Dot Product Constant
    F = math.sqrt(BX * BX + BY * BY + BZ * BZ) * L * math.cos(ang * (math.pi / 180.0))

    ##Constants
    const = math.sqrt(
        math.pow((B * BZ - BY * G), 2)
        * (
            -(F * F) * (A * A + B * B + G * G)
            + (
                B * B * (BX * BX + BZ * BZ)
                + A * A * (BY * BY + BZ * BZ)
                - (2 * A * BX * BZ * G)
                + (BX * BX + BY * BY) * G * G
                - (2 * B * BY) * (A * BX + BZ * G)
            )
            * L
            * L
        )
    )
    denom = (
        (B * B) * (BX * BX + BZ * BZ)
        + (A * A) * (BY * BY + BZ * BZ)
        - (2 * A * BX * BZ * G)
        + (BX * BX + BY * BY) * (G * G)
        - (2 * B * BY) * (A * BX + BZ * G)
    )

    X = (
        (B * B * BX * F) - (A * B * BY * F) + (F * G) * (-A * BZ + BX * G) + const
    ) / denom

    if (B == 0 or BZ == 0) and (BY == 0 or G == 0):
        const1 = math.sqrt(
            G * G * (-A * A * X * X + (B * B + G * G) * (L - X) * (L + X))
        )
        Y = ((-A * B * X) + const1) / (B * B + G * G)
        Z = -(A * G * G * X + B * const1) / (G * (B * B + G * G))
    else:
        Y = (
            (A * A * BY * F) * (B * BZ - BY * G)
            + G * (-F * math.pow(B * BZ - BY * G, 2) + BX * const)
            - A * (B * B * BX * BZ * F - B * BX * BY * F * G + BZ * const)
        ) / ((B * BZ - BY * G) * denom)
        Z = (
            (A * A * BZ * F) * (B * BZ - BY * G)
            + (B * F) * math.pow(B * BZ - BY * G, 2)
            + (A * BX * F * G) * (-B * BZ + BY * G)
            - B * BX * const
            + A * BY * const
        ) / ((B * BZ - BY * G) * denom)

    # Get the new Vector from the origin
    D = Vector(X, Y, Z) + CV
    with warnings.catch_warnings():
        # ignore inconsequential warning
        warnings.simplefilter("ignore")
        temp = calc_dihedral(AV, BV, CV, D) * (180.0 / math.pi)

    di = di - temp
    rot = rotaxis(math.pi * (di / 180.0), CV - BV)
    D = (D - BV).left_multiply(rot) + BV

    return D.get_array()

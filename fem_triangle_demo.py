import argparse
import math


def sub(a, b):
    return [a[0] - b[0], a[1] - b[1]]


def add(a, b):
    return [a[0] + b[0], a[1] + b[1]]


def mul(a, s):
    return [a[0] * s, a[1] * s]


def det2(m):
    return m[0][0] * m[1][1] - m[0][1] * m[1][0]


def inv2(m):
    determinant = det2(m)
    if abs(determinant) < 1e-12:
        raise ValueError("Matrix is singular, cannot invert.")
    inv_det = 1.0 / determinant
    return [
        [m[1][1] * inv_det, -m[0][1] * inv_det],
        [-m[1][0] * inv_det, m[0][0] * inv_det],
    ]


def matmul2(a, b):
    return [
        [
            a[0][0] * b[0][0] + a[0][1] * b[1][0],
            a[0][0] * b[0][1] + a[0][1] * b[1][1],
        ],
        [
            a[1][0] * b[0][0] + a[1][1] * b[1][0],
            a[1][0] * b[0][1] + a[1][1] * b[1][1],
        ],
    ]


def frobenius_norm_sq(m):
    return (
        m[0][0] * m[0][0]
        + m[0][1] * m[0][1]
        + m[1][0] * m[1][0]
        + m[1][1] * m[1][1]
    )


def mat_sub(a, b):
    return [
        [a[0][0] - b[0][0], a[0][1] - b[0][1]],
        [a[1][0] - b[1][0], a[1][1] - b[1][1]],
    ]


def format_vec(v):
    return f"({v[0]:.4f}, {v[1]:.4f})"


def format_mat(m):
    return (
        f"[[{m[0][0]: .4f}, {m[0][1]: .4f}],\n"
        f" [{m[1][0]: .4f}, {m[1][1]: .4f}]]"
    )


def nearest_rotation_2d(F):
    # In 2D, the closest rotation matrix can be recovered from a single angle.
    # This is a compact teaching version of the polar-decomposition idea.
    a = F[0][0] + F[1][1]
    b = F[1][0] - F[0][1]
    theta = math.atan2(b, a)
    c = math.cos(theta)
    s = math.sin(theta)
    return [
        [c, -s],
        [s, c],
    ], theta


class TriangleCase:
    def __init__(self, name, description, X1, X2, X3, x1, x2, x3):
        self.name = name
        self.description = description
        self.X1 = X1
        self.X2 = X2
        self.X3 = X3
        self.x1 = x1
        self.x2 = x2
        self.x3 = x3

    def copied(self):
        return TriangleCase(
            self.name,
            self.description,
            self.X1[:],
            self.X2[:],
            self.X3[:],
            self.x1[:],
            self.x2[:],
            self.x3[:],
        )


def make_case(name, description, X1, X2, X3, x1, x2, x3):
    return TriangleCase(name, description, X1, X2, X3, x1, x2, x3)


CASES = {
    "identity": make_case(
        "identity",
        "No deformation. Current triangle equals reference triangle.",
        [0.0, 0.0], [1.0, 0.0], [0.0, 1.0],
        [0.0, 0.0], [1.0, 0.0], [0.0, 1.0],
    ),
    "stretch_x": make_case(
        "stretch_x",
        "Stretch 2x along x direction while keeping y unchanged.",
        [0.0, 0.0], [1.0, 0.0], [0.0, 1.0],
        [0.0, 0.0], [2.0, 0.0], [0.0, 1.0],
    ),
    "shear": make_case(
        "shear",
        "Shear deformation: the top vertex shifts to the right.",
        [0.0, 0.0], [1.0, 0.0], [0.0, 1.0],
        [0.0, 0.0], [1.0, 0.0], [0.5, 1.0],
    ),
    "compress": make_case(
        "compress",
        "Uniform compression to half size in both directions.",
        [0.0, 0.0], [1.0, 0.0], [0.0, 1.0],
        [0.0, 0.0], [0.5, 0.0], [0.0, 0.5],
    ),
    "rotate": make_case(
        "rotate",
        "Pure 90 degree rotation. This simple demo energy still counts it as deformation.",
        [0.0, 0.0], [1.0, 0.0], [0.0, 1.0],
        [0.0, 0.0], [0.0, 1.0], [-1.0, 0.0],
    ),
}


def edge_matrix(p1, p2, p3):
    e1 = sub(p2, p1)
    e2 = sub(p3, p1)
    return [
        [e1[0], e2[0]],
        [e1[1], e2[1]],
    ]


def triangle_area(p1, p2, p3):
    dm = edge_matrix(p1, p2, p3)
    return abs(det2(dm)) * 0.5


def simple_energy_density(F):
    identity = [[1.0, 0.0], [0.0, 1.0]]
    diff = mat_sub(F, identity)
    return 0.5 * frobenius_norm_sq(diff)


def corotated_energy_density(F):
    R, theta = nearest_rotation_2d(F)
    diff = mat_sub(F, R)
    return 0.5 * frobenius_norm_sq(diff), R, theta


def analyze_case(case):
    Dm = edge_matrix(case.X1, case.X2, case.X3)
    Ds = edge_matrix(case.x1, case.x2, case.x3)
    Dm_inv = inv2(Dm)
    F = matmul2(Ds, Dm_inv)
    area0 = triangle_area(case.X1, case.X2, case.X3)
    psi_simple = simple_energy_density(F)
    psi_corot, R, theta = corotated_energy_density(F)
    energy_simple = area0 * psi_simple
    energy_corot = area0 * psi_corot
    return {
        "Dm": Dm,
        "Ds": Ds,
        "Dm_inv": Dm_inv,
        "F": F,
        "area0": area0,
        "R": R,
        "theta": theta,
        "psi_simple": psi_simple,
        "psi_corot": psi_corot,
        "energy_simple": energy_simple,
        "energy_corot": energy_corot,
    }


def energy_from_case(case, energy_mode):
    result = analyze_case(case)
    if energy_mode == "simple":
        return result["energy_simple"]
    return result["energy_corot"]


def perturb_case(case, node_index, axis_index, delta):
    new_case = case.copied()
    nodes = [new_case.x1, new_case.x2, new_case.x3]
    nodes[node_index][axis_index] += delta
    return new_case


def numerical_force(case, energy_mode, eps=1e-6):
    forces = []
    for node_index in range(3):
        grad = [0.0, 0.0]
        for axis_index in range(2):
            case_plus = perturb_case(case, node_index, axis_index, eps)
            case_minus = perturb_case(case, node_index, axis_index, -eps)
            e_plus = energy_from_case(case_plus, energy_mode)
            e_minus = energy_from_case(case_minus, energy_mode)
            grad[axis_index] = (e_plus - e_minus) / (2.0 * eps)
        forces.append(mul(grad, -1.0))
    return forces


def print_case(case):
    result = analyze_case(case)
    simple_forces = numerical_force(case, "simple")
    corot_forces = numerical_force(case, "corot")

    print(f"Case: {case.name}")
    print(case.description)
    print()
    print("Reference triangle:")
    print(f"  X1 = {format_vec(case.X1)}")
    print(f"  X2 = {format_vec(case.X2)}")
    print(f"  X3 = {format_vec(case.X3)}")
    print()
    print("Current triangle:")
    print(f"  x1 = {format_vec(case.x1)}")
    print(f"  x2 = {format_vec(case.x2)}")
    print(f"  x3 = {format_vec(case.x3)}")
    print()
    print("Dm =")
    print(format_mat(result["Dm"]))
    print()
    print("Ds =")
    print(format_mat(result["Ds"]))
    print()
    print("Dm^-1 =")
    print(format_mat(result["Dm_inv"]))
    print()
    print("F = Ds * Dm^-1 =")
    print(format_mat(result["F"]))
    print()
    print("Nearest rotation R (teaching co-rotational view) =")
    print(format_mat(result["R"]))
    print(f"rotation angle theta = {math.degrees(result['theta']):.2f} degrees")
    print()
    print(f"Reference area A0 = {result['area0']:.4f}")
    print(f"simple psi(F) = 0.5 * ||F - I||_F^2 = {result['psi_simple']:.6f}")
    print(f"simple element energy E = A0 * psi(F) = {result['energy_simple']:.6f}")
    print()
    print(f"co-rotational psi(F) = 0.5 * ||F - R||_F^2 = {result['psi_corot']:.6f}")
    print(f"co-rotational element energy E = A0 * psi(F) = {result['energy_corot']:.6f}")
    print()
    print("Numerical nodal forces from simple energy:")
    print(f"  f1 = {format_vec(simple_forces[0])}")
    print(f"  f2 = {format_vec(simple_forces[1])}")
    print(f"  f3 = {format_vec(simple_forces[2])}")
    print()
    print("Numerical nodal forces from co-rotational energy:")
    print(f"  f1 = {format_vec(corot_forces[0])}")
    print(f"  f2 = {format_vec(corot_forces[1])}")
    print(f"  f3 = {format_vec(corot_forces[2])}")


def print_case_list():
    print("Available FEM triangle cases:")
    for name, case in CASES.items():
        print(f"  {name:10s} - {case.description}")


def main():
    parser = argparse.ArgumentParser(
        description="Minimal 2D triangle FEM teaching demo."
    )
    parser.add_argument(
        "--case",
        choices=sorted(CASES.keys()),
        default="stretch_x",
        help="Which triangle deformation case to analyze.",
    )
    parser.add_argument(
        "--mode",
        choices=["run", "list", "all"],
        default="run",
        help="run analyzes one case, list prints available cases, all analyzes every case.",
    )
    args = parser.parse_args()

    if args.mode == "list":
        print_case_list()
        return

    if args.mode == "all":
        for index, case in enumerate(CASES.values()):
            if index > 0:
                print()
                print("=" * 60)
                print()
            print_case(case)
        return

    print_case(CASES[args.case])


if __name__ == "__main__":
    main()

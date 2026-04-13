import argparse
import math
import os

import matplotlib.pyplot as plt


def add(a, b):
    return [a[0] + b[0], a[1] + b[1]]


def sub(a, b):
    return [a[0] - b[0], a[1] - b[1]]


def mul(a, s):
    return [a[0] * s, a[1] * s]


def length(a):
    return math.sqrt(a[0] * a[0] + a[1] * a[1])


class DistanceConstraint:
    def __init__(self, i, j, rest_length):
        self.i = i
        self.j = j
        self.rest_length = rest_length


class PBDRope:
    def __init__(
        self,
        num_points=8,
        spacing=0.2,
        dt=0.016,
        gravity=(0.0, -9.8),
        num_iterations=10,
        ground_y=-1.5,
        fixed_points=None,
        inverse_masses=None,
    ):
        self.num_points = num_points
        self.spacing = spacing
        self.dt = dt
        self.gravity = [gravity[0], gravity[1]]
        self.num_iterations = num_iterations
        self.ground_y = ground_y
        self.fixed_points = set(fixed_points if fixed_points is not None else [0])
        self.inverse_masses = inverse_masses

        self.positions = []
        self.old_positions = []
        self.velocities = []
        self.inv_masses = []
        self.constraints = []

        self._initialize_particles()
        self._initialize_constraints()

    def _initialize_particles(self):
        for i in range(self.num_points):
            pos = [0.0, -i * self.spacing]
            self.positions.append(pos[:])
            self.old_positions.append(pos[:])
            self.velocities.append([0.0, 0.0])

            if i in self.fixed_points:
                self.inv_masses.append(0.0)
            elif self.inverse_masses is not None:
                self.inv_masses.append(self.inverse_masses[i])
            else:
                self.inv_masses.append(1.0)

    def _initialize_constraints(self):
        for i in range(self.num_points - 1):
            self.constraints.append(DistanceConstraint(i, i + 1, self.spacing))

    def solve_distance_constraint(self, constraint):
        i = constraint.i
        j = constraint.j

        xi = self.positions[i]
        xj = self.positions[j]
        wi = self.inv_masses[i]
        wj = self.inv_masses[j]

        d = sub(xi, xj)
        dist = length(d)
        if dist < 1e-8:
            return

        constraint_error = dist - constraint.rest_length
        direction = mul(d, 1.0 / dist)

        total_inverse_mass = wi + wj
        if total_inverse_mass < 1e-8:
            return

        correction_scale = constraint_error / total_inverse_mass

        correction_i = mul(direction, -wi * correction_scale)
        correction_j = mul(direction, wj * correction_scale)

        self.positions[i] = add(self.positions[i], correction_i)
        self.positions[j] = add(self.positions[j], correction_j)

    def solve_ground_collision(self):
        for i in range(self.num_points):
            if self.positions[i][1] < self.ground_y:
                self.positions[i][1] = self.ground_y

    def step(self):
        for i in range(self.num_points):
            if self.inv_masses[i] == 0.0:
                continue

            self.velocities[i] = add(self.velocities[i], mul(self.gravity, self.dt))
            self.old_positions[i] = self.positions[i][:]
            self.positions[i] = add(self.positions[i], mul(self.velocities[i], self.dt))

        for _ in range(self.num_iterations):
            for constraint in self.constraints:
                self.solve_distance_constraint(constraint)
            self.solve_ground_collision()

        for i in range(self.num_points):
            if self.inv_masses[i] == 0.0:
                continue

            displacement = sub(self.positions[i], self.old_positions[i])
            self.velocities[i] = mul(displacement, 1.0 / self.dt)

    def rope_length_error(self):
        total_error = 0.0
        for constraint in self.constraints:
            i = constraint.i
            j = constraint.j
            dist = length(sub(self.positions[i], self.positions[j]))
            total_error += abs(dist - constraint.rest_length)
        return total_error

    def snapshot(self):
        return {
            "positions": [pos[:] for pos in self.positions],
            "rope_length_error": self.rope_length_error(),
        }

    def print_state(self, frame):
        print(f"frame {frame}")
        for i, pos in enumerate(self.positions):
            print(f"  point {i}: ({pos[0]: .4f}, {pos[1]: .4f})")
        print(f"  total length error: {self.rope_length_error():.6f}")
        print()


def make_config(**overrides):
    config = {
        "name": "default",
        "title": "Default rope",
        "description": "Baseline case: one fixed point, moderate gravity, 10 solver iterations.",
        "num_points": 8,
        "spacing": 0.2,
        "dt": 0.016,
        "gravity": (0.0, -9.8),
        "num_iterations": 10,
        "ground_y": -1.5,
        "fixed_points": [0],
        "inverse_masses": None,
        "num_frames": 120,
    }
    config.update(overrides)
    return config


PRESETS = {
    "default": make_config(),
    "soft": make_config(
        name="soft",
        title="Low iterations: softer rope",
        description="Only 2 solver iterations, so constraints are enforced weakly and the rope appears softer.",
        num_iterations=2,
    ),
    "stiff": make_config(
        name="stiff",
        title="High iterations: stiffer rope",
        description="30 solver iterations, so the rope better preserves segment lengths and looks stiffer.",
        num_iterations=30,
    ),
    "heavy_gravity": make_config(
        name="heavy_gravity",
        title="Stronger gravity",
        description="Gravity is stronger, so the prediction step pulls the rope farther down before projection.",
        gravity=(0.0, -30.0),
    ),
    "two_pins": make_config(
        name="two_pins",
        title="Two fixed points",
        description="Point 0 and point 3 are fixed. The rope shape changes because the constraints must satisfy two anchors.",
        fixed_points=[0, 3],
    ),
    "heavy_tail": make_config(
        name="heavy_tail",
        title="Heavier tail point",
        description="The last point has smaller inverse mass, so it moves less during projection and feels heavier.",
        inverse_masses=[0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.2],
    ),
}


def simulate(config, num_frames=None):
    rope = PBDRope(
        num_points=config["num_points"],
        spacing=config["spacing"],
        dt=config["dt"],
        gravity=config["gravity"],
        num_iterations=config["num_iterations"],
        ground_y=config["ground_y"],
        fixed_points=config["fixed_points"],
        inverse_masses=config["inverse_masses"],
    )

    frames = num_frames if num_frames is not None else config["num_frames"]
    history = []
    for _ in range(frames):
        rope.step()
        history.append(rope.snapshot())
    return rope, history


def print_preset_list():
    print("Available presets:")
    for name, config in PRESETS.items():
        print(f"  {name:12s} - {config['description']}")


def run_text_mode(config, num_frames):
    rope, _ = simulate(config, num_frames=num_frames)
    print("PBD rope demo")
    print(
        f"preset={config['name']}, points={config['num_points']}, spacing={config['spacing']}, "
        f"dt={config['dt']}, gravity={config['gravity']}, iterations={config['num_iterations']}"
    )
    print(f"description: {config['description']}")
    print()

    rope = PBDRope(
        num_points=config["num_points"],
        spacing=config["spacing"],
        dt=config["dt"],
        gravity=config["gravity"],
        num_iterations=config["num_iterations"],
        ground_y=config["ground_y"],
        fixed_points=config["fixed_points"],
        inverse_masses=config["inverse_masses"],
    )

    for frame in range(num_frames):
        rope.step()
        rope.print_state(frame)


def get_bounds(config):
    total_length = (config["num_points"] - 1) * config["spacing"]
    x_padding = 0.3
    y_padding = 0.3
    return (-total_length - x_padding, total_length + x_padding, config["ground_y"] - y_padding, 0.4 + y_padding)


def plot_case(ax, config, history):
    final_positions = history[-1]["positions"]
    initial_positions = [[0.0, -i * config["spacing"]] for i in range(config["num_points"])]

    x0 = [p[0] for p in initial_positions]
    y0 = [p[1] for p in initial_positions]
    x1 = [p[0] for p in final_positions]
    y1 = [p[1] for p in final_positions]

    ax.plot(x0, y0, "o--", color="#b0b0b0", label="initial")
    ax.plot(x1, y1, "o-", color="#d1495b", label="final")
    ax.axhline(config["ground_y"], color="#2b59c3", linestyle=":", label="ground")

    xmin, xmax, ymin, ymax = get_bounds(config)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(config["title"])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True, alpha=0.25)

    final_error = history[-1]["rope_length_error"]
    ax.text(
        0.02,
        0.02,
        f"iterations={config['num_iterations']}\n"
        f"gravity={config['gravity'][1]}\n"
        f"length error={final_error:.5f}",
        transform=ax.transAxes,
        fontsize=9,
        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "#cccccc"},
    )


def save_single_plot(config, history, output_path):
    fig, ax = plt.subplots(figsize=(6, 6))
    plot_case(ax, config, history)
    ax.legend(loc="upper right")
    fig.suptitle(f"PBD Rope Demo: {config['name']}")
    fig.text(0.5, 0.01, config["description"], ha="center", fontsize=10)
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def save_comparison_plot(case_names, output_path):
    cols = 2
    rows = math.ceil(len(case_names) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(12, 5 * rows))
    axes = axes.flatten() if hasattr(axes, "flatten") else [axes]

    for ax, case_name in zip(axes, case_names):
        config = PRESETS[case_name]
        _, history = simulate(config)
        plot_case(ax, config, history)

    for ax in axes[len(case_names):]:
        ax.axis("off")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3)
    fig.suptitle("PBD Rope Parameter Comparison", fontsize=16)
    fig.text(
        0.5,
        0.02,
        "Compare how solver iterations, gravity, fixed points, and inverse mass change the final rope shape.",
        ha="center",
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0.05, 1, 0.93])
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def build_parser():
    parser = argparse.ArgumentParser(
        description="Teaching demo for a minimal PBD rope solver."
    )
    parser.add_argument(
        "--mode",
        choices=["text", "plot", "compare", "list"],
        default="compare",
        help="text prints frames, plot saves one preset figure, compare saves a multi-preset comparison figure, list prints presets.",
    )
    parser.add_argument(
        "--preset",
        choices=sorted(PRESETS.keys()),
        default="default",
        help="Which preset to run in text or plot mode.",
    )
    parser.add_argument(
        "--frames",
        type=int,
        default=None,
        help="Override the number of simulated frames.",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output image path for plot or compare mode.",
    )
    parser.add_argument(
        "--cases",
        nargs="+",
        choices=sorted(PRESETS.keys()),
        default=["default", "soft", "stiff", "heavy_gravity", "two_pins", "heavy_tail"],
        help="Preset names to include in compare mode.",
    )
    return parser


def resolve_output_path(mode, preset, output_path):
    if output_path is not None:
        return output_path

    base_dir = os.path.dirname(os.path.abspath(__file__))
    if mode == "plot":
        return os.path.join(base_dir, f"pbd_{preset}.png")
    if mode == "compare":
        return os.path.join(base_dir, "pbd_compare.png")
    return None


def main():
    parser = build_parser()
    args = parser.parse_args()

    if args.mode == "list":
        print_preset_list()
        return

    config = PRESETS[args.preset]
    num_frames = args.frames if args.frames is not None else config["num_frames"]

    if args.mode == "text":
        run_text_mode(config, num_frames)
        return

    output_path = resolve_output_path(args.mode, args.preset, args.output)
    if args.mode == "plot":
        _, history = simulate(config, num_frames=num_frames)
        save_single_plot(config, history, output_path)
        print(f"Saved plot to: {output_path}")
        print(f"Preset: {config['name']}")
        print(config["description"])
        return

    if args.mode == "compare":
        save_comparison_plot(args.cases, output_path)
        print(f"Saved comparison plot to: {output_path}")
        print("Cases:")
        for case_name in args.cases:
            print(f"  {case_name:12s} - {PRESETS[case_name]['description']}")


if __name__ == "__main__":
    main()

import itertools


def simulate_n_loops(N, loop_range):

    ranges = [
        [1, 2, 3],
        [5, 6, 7],
        [1, 2]
    ]

    print(ranges)

    for combination in itertools.product(*ranges):
        print(combination)


simulate_n_loops(N=3, loop_range=3)
from simulation import biprop_simulate
from cases import full_biprop_case, big_boy_case

from outputs import plot_sim_record, log_results

from matplotlib import pyplot as plt
    

if __name__ == "__main__":

    """
    run case and plot results
    """
    case = big_boy_case

    sim_record = biprop_simulate(case, record=True)

    if sim_record is None:
        raise ValueError("Simulation completed without recording.")
    else:
        print(f"Simulation completed with {len(sim_record.points)} recorded points.")
    
    plot_sim_record(
        simRecord=sim_record,
        file_path=f"plots/{case.name}_test_results.png",
        cols=7,
        show=True
    )

    log_results(
        file_path=f"logs/{case.name}_results.csv",
        sim_record=sim_record
    )
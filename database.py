ID_STATES_PATH = "id_states.json"
SIMULATION_RESULTS_PATH = "simulation_results.csv"


def is_file_exists(path: str) -> bool:
    """ Checks if file with this path exists\n
        input: path to the file\n
        output: true if exists, false otherwise """
    try:
        with open(path, "r") as file:
            return True
    except FileNotFoundError:
        return False


def get_default_id_state() -> dict[str, int]:
    return {"next_gene_id": 0, "next_cog_id": 0, "next_genome_id": 0, "next_ancestor_id" : 0,
            "next_tree_id": 0}


def update_id_state(ids: dict[str,int]):
    """ Updates id_state file with new ids \n
        input: dictionary with all ids needed
        output: none"""
    with open(ID_STATES_PATH, "wt") as file:
        file.write("{\n")
        for i, key in enumerate(ids):
            file.write("\"" + key + "\": " + str(ids[key]))
            if i < len(ids) - 1:
                file.write(",")
            file.write("\n")
        file.write("}")



def get_json_value_from_line(line: str) -> str:
    """ Extracts value from JSON line \n
        input: JSON line \n
        output: value from this line"""
    line = line[line.find(":") + 2: -1]
    if line.find(",") != -1:
        return line[0:line.find(",")]
    return line


def get_id_state() -> dict[str, int]:
    """ Gets ids from id_state file \n
        input: none \n
        output: ids of next gene,COG,genome, ancestor and tree"""
    ids : dict[str, int] = {}
    with open(ID_STATES_PATH, "rt") as file:
        for line in file:
            if line.find(":") == -1:
                continue
            key: str = line[line.find('"') + 1:line.find(":") - 1]
            ids[key] = int(get_json_value_from_line(line))
    return ids


def init_simulation_results():
    """ Initializes simulation results file\n
        input: none \n
        output: none"""
    with open(SIMULATION_RESULTS_PATH, "wt") as file:
        file.write("gene_id,COG_id,genome_id,sequence,ancestor_id,root to leaf distance,tree_id"
                   ",random_seed,ancestor_sequence\n")


def append_simulation_result(result: str, file):
    """ Appends simulation result to simulation results file \n
        input: result\n
        output: none"""
    #with open(SIMULATION_RESULTS_PATH, "a") as file:
    file.write(result + "\n")
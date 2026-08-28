import argparse
import sys

from TreeNode import *
from database import *

id_counter: int = 0


def generate_random_tree(n:int =1) -> TreeNode:
    """ Generates random binary tree structure\n
        input: n - number of leaves\n
        output: root of the tree """
    node:TreeNode = TreeNode(Genome())
    if n == 1:
        global id_counter
        node.leaf_id =id_counter
        id_counter+= 1
        return node
    leftLeaves:int = random.randint(1,n - 1)
    rightLeaves:int = n - leftLeaves
    node.left = generate_random_tree(leftLeaves)
    node.right = generate_random_tree(rightLeaves)
    return node


def generate_edge_lengths_exponential(current: TreeNode|None = None, mean: float= 0.1):
    """ Generates lengths of edges using exponential distribution
        input: current tree node and mean edge length
        output: none"""
    if current is None:
        return
    _lambda: float = 1.0 / mean
    while True:
        uniform: float = random.random()
        current.edge_length = -math.log(1 - uniform) / _lambda
        if EvolutionTree.get_probability_from_length_jukes_cantor(current.edge_length) < 0.5:
            break
    generate_edge_lengths_exponential(current.left, mean)
    generate_edge_lengths_exponential(current.right, mean)


def generate_simulation_seed() -> int:
    """Generates random seed\n
        input: none\n
        output: random seed """
    return random.randint(0,sys.maxsize)


def summarise_simulation(cog_id: int, ancestor_id: int, tree_id: int, seed: str,
                         current: TreeNode, distance: float,ids: dict[str, int], root:TreeNode,
                         file):
    """writes results of the simulation into file\n
        input: ids of COG, ancestor, tree, seed that was used for simulation,
        current node in the tree, distance from the root next ids root node and results file\n
        output: none
        """
    if current is None:
        return
    distance += current.edge_length
    if current.left is None and current.right is None: #leaf
        genome_id = ids["next_genome_id"]
        ids["next_genome_id"] = genome_id + 1
        for i, gene in enumerate(current.val.genes):
            gene_id = ids["next_gene_id"]
            ids["next_gene_id"] = gene_id + 1
            append_simulation_result("gene_" + str(gene_id) + ",COG_" + str(cog_id + i) + ",G_" +
                                     str(genome_id) + "," + gene + ",ancestor_" + str(ancestor_id + i)
                                     + "," + str(distance) + ",tree_" + str(tree_id) + "," + seed +
                                     "," + root.val.genes[i],file)
    summarise_simulation(cog_id, ancestor_id, tree_id, seed, current.left, distance, ids, root, file)
    summarise_simulation(cog_id, ancestor_id, tree_id, seed, current.right, distance, ids, root, file)


def get_arguments():
    parser = argparse.ArgumentParser(description="Generate random mutation tree")
    parser.add_argument("--n", help="number of leaves", default=1000, type=int)
    parser.add_argument("--m", help="number of COGs", default=200, type=int)
    parser.add_argument("--l", help="length of gene", default=10, type=int)
    parser.add_argument("--mean", help="mean edge length", default=0.1, type=float)
    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    #generating structure of a tree and its edges
    root:TreeNode = generate_random_tree(args.n)
    tree: EvolutionTree = EvolutionTree(root, args.n)
    print("finished generating structure")

    generate_edge_lengths_exponential(root.left,args.mean)
    generate_edge_lengths_exponential(root.right, args.mean)
    print("finished generating lengths")

    tree.calculate_leaf_distances()

    #setting random seed
    random_seed = generate_simulation_seed()
    random.seed(random_seed)

    if not is_file_exists(ID_STATES_PATH):
        update_id_state(get_default_id_state())
    if not is_file_exists(SIMULATION_RESULTS_PATH):
        init_simulation_results()
    #setting tree id
    ids: dict[str,int] = get_id_state()
    tree.tree_id = ids["next_tree_id"]
    ids["next_tree_id"] = tree.tree_id + 1

    #beginning of simulation
    root.val.generate_random_genome(args.m,args.l)
    tree.ancestor_id = ids["next_ancestor_id"]
    ids["next_ancestor_id"] = tree.ancestor_id + args.m
    tree.cog_id = ids["next_cog_id"]
    ids["next_cog_id"] = tree.cog_id + args.m

    tree.mutate_tree_jukes_cantor()
    print("finished mutating tree")
    with open(SIMULATION_RESULTS_PATH, "a") as file:
        summarise_simulation(tree.cog_id, tree.ancestor_id, tree.tree_id, str(random_seed),
                         root,0,ids,tree.root, file)
    print("finished writing results in table")
    #print(root.convert_to_newick_format(True) + ";")
    #tree.print_tree_dfs()
    tree.write_tree_newick_format("tree_" + str(tree.tree_id) + ".tree")

    update_id_state(ids)


if __name__ == "__main__":
    main()
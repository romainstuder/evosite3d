#!/usr/bin/env python3
"""Parse the NCBI taxonomy"""

import logging
import random
from typing import Dict, List, Tuple

LOGGER = logging.getLogger()
LOGGER.setLevel(0)


class Node:
    """Definition of the class Node"""

    def __init__(self):
        self.tax_id = "0"  # Number of the tax id.
        self.parent = "0"  # Number of the parent of this node
        self.children = []  # List of the children of this node
        self.division = None  # Division.
        self.is_tip = True  # Tip = True if it's a terminal node, False if not.
        self.name = ""  # Name of the node: taxa if it's a terminal node, number if not.


def get_genealogy(name_object, leaf_node: str) -> List[str]:
    """Trace genealogy from root to leaf"""
    ancestors = []  # Initialise the list of all nodes from root to leaf.
    gen_tax_id = leaf_node  # Define leaf
    while 1:
        if gen_tax_id in name_object:
            ancestors.append(gen_tax_id)
            gen_tax_id = name_object[gen_tax_id].parent  # Move up to parents
        else:
            break
        if gen_tax_id == "1":
            # If it is the root, we reached the end.
            # Add it to the list and break the loop
            ancestors.append(gen_tax_id)
            break
    return ancestors  # Return the list


def _get_all_descendant_nodes(name_object, taxid: str) -> List[str]:
    """Get all descendant of a node"""
    descendant_nodes: List[str] = [taxid]
    if len(name_object[taxid].children) > 0:
        for child in name_object[taxid].children:
            descendant_nodes = descendant_nodes + _get_all_descendant_nodes(name_object, child)
    return descendant_nodes


def _keep_terminal(name_object, nodes_list) -> List[str]:
    """Keep only terminal nodes"""
    return [x for x in nodes_list if name_object[x].is_tip]


def _keep_division(name_object, nodes_list, target_division) -> List[str]:
    """Keep only division nodes"""
    return [x for x in nodes_list if name_object[x].division == target_division]


def get_all_descendants(name_object, target_division: str, taxid: str) -> List[str]:
    """Get all taxa of a node"""

    terminal_nodes = _get_all_descendant_nodes(name_object, taxid)
    terminal_nodes = _keep_division(name_object, terminal_nodes, target_division)

    return terminal_nodes  # Return a list


def get_common_ancestor(name_object, node_list: List[str]):
    """Function to find common ancestor between two nodes or more

    Args:
        name_object (name_object): taxonomy to use
        node_list (list): list of node

    Returns:
        node (str): node of the common ancestor between nodes
    """

    # global name_object
    list1 = get_genealogy(name_object, node_list[0])  # Define the whole genealogy of the first node
    ancestral_list: List[str] = []
    for node in node_list:
        list2 = get_genealogy(name_object, node)  # Define the whole genealogy of the second node
        ancestral_list = []
        for taxid in list1:
            if taxid in list2:  # Identify common nodes between the two genealogy
                ancestral_list.append(taxid)
        list1 = ancestral_list  # Reassigning ancestral_list to list 1.
    last_common_ancestor = ancestral_list[
        0
    ]  # Finally, the first node of the ancestral_list is the common ancestor of all nodes.
    return last_common_ancestor  # Return a node


def load_ncbi_names(filename: str = "names.dmp") -> Tuple[Dict, Dict]:
    """Load NCBI names definition ("names.dmp")

    Args:
        filename (str): filename of NCBI names

    Returns:
        name_dict, name_dict_reverse

    """

    name_dict = {}  # Initialise dictionary with TAX_ID:NAME
    name_dict_reverse = {}  # Initialise dictionary with NAME:TAX_ID

    LOGGER.warning(f"Load {filename}")
    name_file = open(filename, "r")
    while 1:
        line = name_file.readline()
        if line == "":
            break
        line = line.rstrip()
        line = line.replace("\t", "")
        tab = line.split("|")
        if tab[3] == "scientific name":
            tax_id, name = tab[0], tab[1]  # Assign tax_id and name ...
            name_dict[tax_id] = name  # ... and load them
            name_dict_reverse[name] = str(tax_id)  # ... into dictionaries
    name_file.close()
    return name_dict, name_dict_reverse


def load_ncbi_taxonomy(name_dict, filename: str = "nodes.dmp"):
    """Load taxonomy NCBI file ("nodes.dmp")

    Args:
        filename (str): filename of ncbi taxonomy
        name_dict (dict): name_dict

    Returns:

    """

    # Define taxonomy variable
    # global name_object
    name_object: Dict = {}

    LOGGER.warning(f"Load {filename}")
    taxonomy_file = open(filename, "r")
    while 1:
        line = taxonomy_file.readline()
        if line == "":
            break
        line = line.replace("\t", "")
        tab = line.split("|")

        tax_id = str(tab[0])
        tax_id_parent = str(tab[1])
        division = str(tab[2])

        # Define name of the taxonomy id
        name = "unknown"
        if tax_id in name_dict:
            name = name_dict[tax_id]

        if tax_id not in name_object:
            name_object[tax_id] = Node()
        name_object[tax_id].tax_id = tax_id  # Assign tax_id
        name_object[tax_id].parent = tax_id_parent  # Assign tax_id parent
        name_object[tax_id].name = name  # Assign name
        name_object[tax_id].division = division  # Assign name

        # Add it has children to parents
        children_list = []
        if tax_id_parent in name_object:
            children_list = name_object[tax_id_parent].children  # If parent is in the object
        else:
            name_object[tax_id_parent] = Node()
            name_object[tax_id_parent].tax_id = tax_id_parent  # Assign tax_id
        children_list.append(tax_id)  # ... we found its children.
        name_object[tax_id_parent].children = children_list  # ... so add them to the parent

        # As the parent node is found, it is not a terminal node then
        name_object[tax_id_parent].is_tip = False

    taxonomy_file.close()

    return name_object


def main():
    """Main function"""

    # Load name_dict, name_dict_reverse and taxonomy
    name_dict, name_dict_reverse = load_ncbi_names(filename="names.dmp")  # Load names
    ncbi_taxonomy = load_ncbi_taxonomy(filename="nodes.dmp", name_dict=name_dict)

    #################################################
    #                                               #
    #    Example 1 : Evolutionary history of human  #
    #                                               #
    #################################################

    LOGGER.warning('Load taxonomy NCBI file ("nodes.dmp")')

    # But what is the tax_id of Human ???
    tax_id_human = name_dict_reverse["Homo sapiens"]

    LOGGER.warning(f"Tax_id of Human (Homo sapiens) is: {tax_id_human}")

    # Ok, now define human genealogy...
    human_genealogy = get_genealogy(ncbi_taxonomy, tax_id_human)
    print("human_genealogy", human_genealogy)
    # ... and display it, with tax_id and name
    for tax_id in human_genealogy:
        print(
            f"Name: {ncbi_taxonomy[tax_id].name}, "
            f"Tax_id: {tax_id}, "
            f"Division: {ncbi_taxonomy[tax_id].division}"
        )

    ##################################################################
    #                                                                #
    #    Example 2 : Common ancestor  between Trichoplax and human   #
    #                                                                #
    ##################################################################

    # What is the common ancestor between Trichoplax and Human?

    # Trichoplax: http://en.wikipedia.org/wiki/Trichoplax

    # Define the two nodes and add them to a list
    tax_id_1 = name_dict_reverse["Trichoplax adhaerens"]
    tax_id_2 = name_dict_reverse["Homo sapiens"]
    list_of_nodes = [tax_id_1, tax_id_2]
    print("list_of_nodes", list_of_nodes)

    # Identify the common ancestor
    common_ancestor = get_common_ancestor(ncbi_taxonomy, list_of_nodes)

    print(
        f"The common ancestor between {ncbi_taxonomy[tax_id_1].name} and "
        f"{ncbi_taxonomy[tax_id_2].name} is: {ncbi_taxonomy[common_ancestor].name}"
    )

    ###################################################################
    #                                                                 #
    #    Example 3 : Get all species (terminal nodes) of a node       #
    #                                                                 #
    ###################################################################

    name = "Primates"
    tax_id = name_dict_reverse[name]
    target_division = "genus"
    print(name, tax_id)

    all_descendants = _get_all_descendant_nodes(ncbi_taxonomy, tax_id)
    all_descendants = _keep_division(ncbi_taxonomy, all_descendants, target_division)

    print(f"Number of all_species of {name}: {len(all_descendants)}")

    sample_list = all_descendants
    if len(all_descendants) > 20:
        sample_list = random.sample(all_descendants, 20)

    for tax_id in sample_list:
        print(
            f"Name: {ncbi_taxonomy[tax_id].name}, "
            f"Tax_id: {tax_id}, "
            f"Division: {ncbi_taxonomy[tax_id].division}, "
            f"Is tip: {ncbi_taxonomy[tax_id].is_tip}, "
            f"Childrien: {ncbi_taxonomy[tax_id].children} "
        )


if __name__ == "__main__":
    main()

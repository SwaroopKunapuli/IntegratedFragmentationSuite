import itertools
import os

import pypackmol as pyp

class Moiety:
    def __init__(self, name: str, charge: float, xyz_file: str) -> None:
        self.name = name
        self.charge = charge
        self.xyz_file = xyz_file
    
    def __str__(self) -> str:
        return(f"Name: {self.name}\t Charge: {self.charge}\t XYZ: {self.xyz_file}")
        
class Fragment:
    def __init__(self, id: list, fundamental_moieties:list, graph_path=str, level:str = None, charge: float = None, parents: list = None, daughters: list = None) -> None:
        self.id = id
        self.moieties = fundamental_moieties
        if len(self.id) != len(self.moieties):
            raise ValueError("The number of moieties and the id do not match")
        for i in range(len(self.id)):
            if self.id[i][1] != self.moieties[i].name:
                raise ValueError(f"The moiety {self.id[i][1]} does not match the moiety {self.moieties[i].name}")
        self.level = level if level != None else sum([x[0] for x in self.id])
        self.charge = charge if charge != None else self.compute_charge()
        self.parents = parents if parents != None else self.get_fragment_parents()
        self.daughters = daughters if daughters != None else self.get_fragment_daughters()
        self.data_location = f"{graph_path}/Lvl_{self.level}/{format_id(self.id, sep='_')}"
        self.prepare_directory()
        #self.gen_xyz()
        
    def __str__(self) -> str:
        parents = '\n'.join([display_id(parent) for parent in self.parents])
        daughters = '\n'.join([display_id(daughter) for daughter in self.daughters])
        return(f"###########################\nLEVEL\t{self.level}\nCHARGE\t{self.charge}\nID\n{display_id(self.id)}\nPARENTS\n{parents}\nDAUGHTERS\n{daughters}\n")
    
    def __iter__(self) -> list:
        return iter(self.id)
    
    def compute_charge(self) -> float:
        charge = 0
        for i in range(len(self.id)):
            charge += self.id[i][0] * self.moieties[i].charge
        return(charge)
    
    def get_fragment_parents(self) -> list:
        parents = []
        for i in range(len(self.id)):
            temp_id = self.id.copy()
            temp_id[i] = (temp_id[i][0]+1, temp_id[i][1])
            parents += [temp_id]
        return(parents)
    
    def get_fragment_daughters(self) -> list:
        daughters = []
        for i in range(len(self.id)):
            if self.id[i][0] > 0:
                temp_id = self.id.copy()
                temp_id[i] = (temp_id[i][0]-1, temp_id[i][1])
                daughters += [temp_id]
        return(daughters)
    
    def generate_xyz(self):
        pm = pyp.Packmol(dimension=10)
        for fragment_composition, moiety_object in zip(self, self.moieties):
            if fragment_composition[1] == moiety_object.name:
                if fragment_composition[0] > 0:
                    pm.add_structure(moiety_object.xyz_file,count=fragment_composition[0],input_format="xyz")
            else:
                raise ValueError(f"Moiety order mismatch for fragment {self.id}")
        result = pm.pack(output=f"{self.data_location}/{format_id(self.id, sep='_')}.xyz")

    def prepare_directory(self) -> None:
        if os.path.exists(self.data_location):
            pass
        else:
            os.makedirs(self.data_location)
    
    def prepare_crest(self) -> None:
        os.chdir(self.data_location)
        os.makedirs("CREST")
        os.chdir("CREST")
    
    def run_crest(self) -> None:
        pass
    
    def process_crest(self) -> None:
        pass
        
class TopParent(Fragment):
    def __init__(self, composition: list, level:str = None, charge: float = None, parents: list = None, daughters: list = None) -> None:
        id = [(composition[i][0], composition[i][1].name) for i in range(len(composition))]
        fundamental_moieties = [moiety[1] for moiety in composition]
        self.top_parent = True
        super().__init__(id, fundamental_moieties, level, charge, parents, daughters)
        self.level = sum([x[0] for x in self.id])
        self.parents = []
        
    def __str__(self) -> str:
        daughters = '\n'.join([display_id(daughter) for daughter in self.daughters])
        return(f"###########################\nLEVEL\t{self.level}\nCHARGE\t{self.charge}\nID\n{display_id(self.id)}\nPARENTS:\nThis is the top parent\nDAUGHTERS\n{daughters}\n")
    
    def prepare_directory(self) -> None:
        return None
    
    def gen_xyz(self) -> None:
        return None

class Fragment_TopParentGraph(Fragment):
    ### We need to exclude parents that are not compatible with the top parent ion composition when creating the parent list
    def __init__(self, id: list, fundamental_moieties:list, graph_path:str, top_parent:TopParent, level:str = -1, charge: float = None, parents: list = None, daughters: list = None) -> None:
        self.top_parent = top_parent
        level = sum([moiety[0] for moiety in id])
        super().__init__(id, fundamental_moieties, graph_path, level, charge, parents, daughters)
        
    def get_fragment_parents(self) -> list:
        parents = []
        for i in range(len(self.id)):
            temp_id = self.id.copy()
            temp_id[i] = (temp_id[i][0]+1, temp_id[i][1])
            if temp_id[i][0] > self.top_parent.id[i][0]:
                pass
            else:
                parents += [temp_id]
        return(parents)
    
class FragmentationGraph:
    def __init__(self, fundamental_moieties: list, graph_path: str=os.getcwd()) -> None:
        self.moieties = fundamental_moieties
        self.id_pattern = [(0,moiety.name) for moiety in fundamental_moieties]
        self.data_path = graph_path
        self.lvl_fragments = []
    
    def plot(self) -> None:
        try: 
            import networkx as nx
            import matplotlib.pyplot as plt
        except:
            print("Please install networkx and matplotlib to plot the fragmentation tree")
            return
        edges = []
        G = nx.DiGraph()
        for idx_lvl, lvl in enumerate(self.lvl_fragments):
            for fragment in lvl:
                G.add_node(format_id(fragment.id), layer=idx_lvl)
                if idx_lvl < len(self.lvl_fragments)-1:
                    for parent in fragment.parents:
                        edges += [(format_id(fragment.id), format_id(parent))]
                    G.add_edges_from(edges)
        pos = nx.multipartite_layout(G, subset_key='layer', align='horizontal', scale=100)
        nx.draw(G, pos, with_labels=True, node_shape='s', node_size=1, font_weight='bold', font_size=6, bbox = dict(facecolor = "skyblue"))
        plt.show()   
    
class FragmentationGraph_Full(FragmentationGraph):
    def __init__(self, fundamental_moieties: list, max_size: int, graph_path: str = os.getcwd()) -> None:
        super().__init__(fundamental_moieties, graph_path)
        self.lvl_fragments = self.generate_fragments(max_size)
        self.size = len(self.lvl_fragments)
    
    def __str__(self) -> str:
        formatted_moieties = '\t'.join([str(moiety.name) for moiety in self.moieties])
        return(f"Levels: {self.size}\nComposition: {formatted_moieties}\nTop level fragments count: {len(self.lvl_fragments[-1])}\n")
        
    def generate_fragments(self, size: int) -> None:
        lvl_fragments = []
        moiety_names = [moiety.name for moiety in self.moieties]
        for lvl in range(1, size+1):
            combos = [x for x in itertools.combinations_with_replacement(moiety_names, lvl)]
            fragments = []
            for el in combos:
                id = get_fragment_ID(el, self.id_pattern)
                fragments += [Fragment(id, self.moieties, self.data_path, lvl)]
            lvl_fragments += [fragments]
        return(lvl_fragments)
        
class FragmentationGraph_TopParent(FragmentationGraph):
    def __init__(self, parent_ion: TopParent, graph_path = os.getcwd()) -> None:
        fundamental_moieties = parent_ion.moieties
        self.top_parent = parent_ion
        super().__init__(fundamental_moieties, graph_path)
        self.lvl_fragments = self.generate_fragments()
        
    def __str__(self) -> str:
        moieties = '\t'.join([str(moiety.name) for moiety in self.moieties])
        return(f"Levels: {self.top_parent.level}\nComposition: {moieties}\nTop parent ion: {self.top_parent.id}\n")
        
    def generate_fragments(self) -> None:
        lvl_fragments = []
        top_parent = self.top_parent
        top_lvl = 0
        for moiety in top_parent:
            top_lvl += moiety[0]
        moiety_names = [moiety.name for moiety in self.moieties]
        for n in range(1, top_lvl+1):
            combos = [x for x in itertools.combinations_with_replacement(moiety_names, n)]
            fragments = []
            for el in combos:
                for moiety in top_parent:
                    if el.count(moiety[1]) > moiety[0]:
                        break
                else:
                    id = get_fragment_ID(el, self.id_pattern)
                    fragments += [Fragment_TopParentGraph(id, self.moieties, self.data_path, top_parent)]
            lvl_fragments += [fragments]
        return(lvl_fragments)

### Global functions ###
        
def get_fragment_ID(fundamental_moieties: list, id_pattern: list = None) -> list:
    id = []
    if id_pattern != None:
        for i in range(len(id_pattern)):
            if id_pattern[i][1] in fundamental_moieties:
                id += [(fundamental_moieties.count(id_pattern[i][1]),id_pattern[i][1])]
            else:
                id += [(0,id_pattern[i][1])]
    else:
        for moiety in set(fundamental_moieties):
            id += [(fundamental_moieties.count(moiety),moiety)]
    return(id)

def display_id(id: list) -> str:
    return("".join([f"{x[0]} {x[1]} + " for x in id])[:-3])

def format_id(id: list, sep=" ") -> str:
    return(sep.join([f"{x[0]}{x[1][0:2]}" for x in id])[:])

###### TESTING FUNCTIONS ######
def test_graph_full(moiety_list: list, max_size: int) -> None:
    graph_full = FragmentationGraph_Full(moiety_list, max_size, "./test_full")
    graph_full.plot()
    
def test_graph_top_parent(moiety_list: list, parent_ion: list) -> None:
    graph_with_top_parent = FragmentationGraph_TopParent(parent_ion, "./test_parent")
    graph_with_top_parent.plot()
    
if __name__ == "__main__":
    # Building the moieties objects needed about everywhere
    Cl = Moiety("Cl", -1, "/home/simon/Bureau/Softwares/Fragmentation_Code/Cl.xyz")
    Urea = Moiety("Urea", 0, "/home/simon/Bureau/Softwares/Fragmentation_Code/Urea.xyz")
    Choline = Moiety("Choline", 1, "/home/simon/Bureau/Softwares/Fragmentation_Code/Choline.xyz")
    moiety_list = [Cl, Urea, Choline]
    
    ## Testing things
    test_graph_full(moiety_list, 10)
    top_parent_ion = TopParent(composition=[(6, Cl), (5, Urea), (4, Choline)])
    test_graph_top_parent(moiety_list, top_parent_ion)

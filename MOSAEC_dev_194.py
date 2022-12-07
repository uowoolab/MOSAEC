import os
import pandas as pd
from ccdc import io
from ccdc import crystal
from ccdc import molecule
import mendeleev
from ccdc import conformer
from ccdc import descriptors
import math

def readentry (input_cif):
    """reads the input cif into a CCDC.Crystal type object,
    and converts to standard atom labelling convention
    (necessary for the rest of the code to function properly"""
    
    #read in the cif to a crystal object
    with io.CrystalReader(input_cif, format='cif') as readcif:
        cif = readcif[0]
    readcif.close()

    # to remove duplicate atoms, need the empirical formula
    formula = cif.formula
    elamnt = formula.split(' ')

    #now convert to standard labelling convention and identify
    #duplicate atoms to be removed
    with open(input_cif, 'r') as file:
        file.seek(0)
        newstring = str()
        lines = file.readlines()
        loop_pos = 0
        start = 0
        end = 0
        columncount = 0
        type_pos = 0
        label_pos = 0
        x_pos = 0
        y_pos = 0
        z_pos = 0
        for i, line in enumerate(lines):
            lines[i] = lines[i].lstrip()
        for i, line in enumerate(lines):
            #locate atom type and site label columns
            if 'loop_' in line:
                loop_pos = i
            if ("_atom" in line) and (not "_geom" in line) and (not "_aniso" in line):
                start = loop_pos + 1
                end = i + 1
        for i in range(start, end):
            if 'atom_site_type_symbol' in lines[i]:
                type_pos = columncount
            if 'atom_site_label' in lines[i]:
                label_pos = columncount
            columncount += 1
        counting = {}
        cutoff = {}
        to_remove = []
        for i in range(end, len(lines)):
            if 'loop_' in lines[i]:
                break
            #lines with atom information will contain a ., so only look at these
            if '.' in lines[i]:
                # split lines by whitespace
                col = lines[i].split()
                #keep count of how many of each element type
                if not col[type_pos] in counting:
                    counting[col[type_pos]] = 1
                elif col[type_pos] in counting:
                    counting[col[type_pos]] +=1
                #new atom labels
                newlabel = f'{col[type_pos]}{counting[col[type_pos]]}'
                lines[i] = lines[i].replace(col[label_pos], newlabel)
                #cutoff repeated atoms
                if newlabel in elamnt:
                    cutoff[col[type_pos]] = counting[col[type_pos]]                
                if col[type_pos] in cutoff:
                    if counting[col[type_pos]] > cutoff[col[type_pos]]:
                        to_remove.append(lines[i])
        #remove unnecessary atoms
        for i in to_remove:
            lines.remove(i)
        #combine to new string
        for i in lines:
            newstring += i
        #read into new crystal object and assign bonds
        newcif = crystal.Crystal.from_string(newstring, format='cif')
        newcif.assign_bonds()
        file.close()   
    return (newcif)
    
def readSBU (input_mol2):
    """MOSAEC cnan now handle SBUs, but expectes them to be in .mol2"""
    #First have to convert the connection points (unknown atom types)
    #to H. Has to be done as a string so CSD reader can 
    #interpret bonding
    with open(input_mol2, 'r') as file:
        file.seek(0)
        newstring = str()
        lines = file.readlines()
        for i, line in enumerate(lines):
            if '*' in line:
                lines[i] = lines[i].replace('Du', 'H')
            newstring += lines[i]    
        cif = crystal.Crystal.from_string(newstring, format='mol2')
        cif.assign_bonds()
        file.close()

    mol = cif.molecule
    #MOSAEC needs atoms to have unique labels, so make them unique
    count = 1
    for atom in mol.atoms:
        atom.label = f'{atom.label}{count}'
        count += 1

    return(mol)

def read_CSD_entry (input_cif):
    """read entries directly from the CSD"""
    #read in the cif to a crystal object
    csd_crystal_reader = io.CrystalReader('CSD')
    cif = csd_crystal_reader.crystal(input_cif)
    cif.assign_bonds()
    csd_crystal_reader.close()
    return (cif)

def get_no_metal_molecule(inputmolecule):
    workingmol = inputmolecule.copy()
    for atom in workingmol.atoms:
        if atom.is_metal:
            workingmol.remove_atom(atom)
    workingmol.assign_bond_types(which = 'All')
    return(workingmol)

def get_unique_sites(mole, asymmole):
    #blank list for unique sites
    uniquesites = []
    labels = []
    asymmcoords = []
    molecoords = []
    duplicates = []
    for atom in asymmole.atoms:
        asymmcoords.append(atom.coordinates)
    for atom in mole.atoms:
        if atom.coordinates in asymmcoords:
            if not atom.coordinates in molecoords:
                if not atom.label in labels:
                    uniquesites.append(atom)
                    molecoords.append(atom.coordinates)
                    labels.append(atom.label)
                else:
                    duplicates.append(atom)
            else:
                duplicates.append(atom)
    if len(duplicates) >= 1:
        for datom in duplicates:
            for atom in uniquesites:
                if any([
                    (datom.coordinates == atom.coordinates),
                    (datom.label == atom.label)
                    ]):
                    if datom.atomic_symbol == atom.atomic_symbol:
                        if len(datom.neighbours) > len(atom.neighbours):
                            uniquesites.remove(atom)
                            uniquesites.append(datom)
                    if not datom.label in labels:
                        uniquesites.append(datom)
                        labels.append(datom.label)
    return (uniquesites)

def get_metal_sites(sites):
    metalsites = []
    for site in sites:
        if site.is_metal == True:
            metalsites.append(site)
    return metalsites

def get_ligand_sites(metalsites,sites):
    metal_sphere = {}
    for metal in metalsites:
        sphere1 = []
        for ligand in metal.neighbours:
            if not ligand.is_metal == True:
                for site in sites:
                    if ligand.label == site.label:
                        sphere1.append(site)
        metal_sphere[metal] = sphere1
    return metal_sphere

def get_binding_sites(metalsites, uniquesites):
    binding_sites = set()
    for metal in metalsites:
        for ligand in metal.neighbours:
            for site in uniquesites:
                if ligand.label == site.label:
                    binding_sites.add(site)
    return (binding_sites)

def get_sphere_2(sphere_1, sites):
    sphere2 = {}
    for metal in sphere_1:
        for ligand in sphere_1[metal]:
            sphere = []
            for environ in ligand.neighbours:
                for site in sites:
                    if environ.label == site.label:
                        sphere.append(site)
            sphere2[ligand] = sphere
    return sphere2

def ringVBOs(mole):
    ringVBO = {}
    unassigned = mole.atoms
    ringcopy = mole.copy()
    oncycle_atoms = []
    offcycle_atoms = []
    oncycle_labels = []
    offcycle_labels = []

    #remove all the metals, this 
    #prevents metal-containing rings (i.e. pores)
    #from interfering
    for atom in ringcopy.atoms:
        if atom.is_metal:
            ringcopy.remove_atom(atom)
    
    #collect all the cyclic atoms
    for atom in ringcopy.atoms:
            if atom.is_cyclic:
                if not atom in oncycle_atoms:
                    oncycle_atoms.append(atom)
                    oncycle_labels.append(atom.label)
    
    #we also need everything that the cyclic atoms are bound to
    for atom in oncycle_atoms:
        for neighbour in atom.neighbours:
            if not neighbour in oncycle_atoms:
                if not neighbour in offcycle_atoms:
                    offcycle_atoms.append(neighbour)
                    offcycle_labels.append(neighbour.label)

    cyclicsystem = (oncycle_atoms + offcycle_atoms)
    
    #remove every atom that isn't part of or directly bound to a cycle
    for atom in ringcopy.atoms:
        if not atom in cyclicsystem:
            ringcopy.remove_atom(atom)
    
    #find all non-cyclic bonds
    #bonds between cycles, break and cap with H
    for bond in ringcopy.bonds:
        if not bond.is_cyclic:
            #bonds between cycles
            if all((member.label in oncycle_labels for member in bond.atoms)):
                member1 = bond.atoms[0]
                member2 = bond.atoms[1]
                Hcap1 = molecule.Atom('H', coordinates = member1.coordinates)
                Hcap2 = molecule.Atom('H', coordinates = member2.coordinates)
                Hcap1_id = ringcopy.add_atom(Hcap1)
                Hcap2_id = ringcopy.add_atom(Hcap2)
                ringcopy.add_bond(bond.bond_type, Hcap1_id, member2)
                ringcopy.add_bond(bond.bond_type, Hcap2_id, member1)
                ringcopy.remove_bond(bond)
    
    #cap off-cycle atoms
    for offatom in offcycle_atoms:
        for bond in offatom.bonds:
            if bond.bond_type == 'Single':
                offatom.atomic_symbol = 'H'
            elif bond.bond_type == 'Double':
                offatom.atomic_symbol = 'O'
                
    #for each cyclic system, reassign bonds, kekulize, and get VBO
    #the bond and atom pruning we did above ensures that fused cycles
    #will be treated as a single system
    #while non-fused cycles that are connected via bonding are treated
    #as seperate systems
    for cyclesys in ringcopy.components:
        #reassign bonds and kekulize
        cyclesys.assign_bond_types()
        cyclesys.kekulize()

        #quick check for delocalized systems in the ring
        #if there are any, get the delocalised bond orders
        if any (bond.bond_type == 'Delocalised' for bond in cyclesys.bonds):
            rdVBO = delocalisedLBO(cyclesys)

        #assign VBO for each on-cycle atom
        for ratom in cyclesys.atoms:
            rVBO = 0
            if ratom.label in oncycle_labels:
                for rbond in ratom.bonds:
                #Each bond contributes to Ligand Bond Order according to its type
                    if rbond.bond_type == 'Single':
                        rVBO +=1
                    elif rbond.bond_type == 'Double':
                        rVBO +=2
                    elif rbond.bond_type == 'Triple':
                        rVBO +=3
                    elif rbond.bond_type == 'Quadruple':
                        rVBO +=4
                    elif rbond.bond_type == 'Delocalised':
                        rVBO += rdVBO[ratom]
                    elif rbond.bond_type == 'Aromatic':
                        rVBO += 0
                        print ("impossible Aromatic bond")
                
                #the VBOs are currently associated to atom objects
                #in molecule objects that we have modified
                #we need these to be associated to atom objects in
                #the parent (unmodified) molecule object
                for matom in unassigned:
                    if matom.label == ratom.label:
                        ringVBO[matom] = rVBO
                        unassigned.remove(matom)
    return(ringVBO)


def assign_VBS(atom, rVBO, dVBO):
    """This function will assign a Valence-Bond-Sum (VBS) to an atom. 
    Takes one CCDC atom object, list of binding sites, and the metal-free
    molecule object as inputs"""
    VBO = 0
    if atom.is_metal:
        return(0)
    if atom in rVBO:
        VBO = rVBO[atom]
    else:
        for bond in atom.bonds:
            if any(batom.is_metal
            for batom in bond.atoms):
                VBO += 0
            #Each bond contributes to Ligand Bond Order according to its type
            elif bond.bond_type == 'Single':
                VBO +=1
            elif bond.bond_type == 'Double':
                VBO +=2
            elif bond.bond_type == 'Triple':
                VBO +=3
            elif bond.bond_type == 'Quadruple':
                VBO +=4
            elif bond.bond_type == 'Delocalised':
                VBO += dVBO[atom]
            elif bond.bond_type == 'Aromatic':
                VBO += aVBO[atom]
    return(VBO)

def delocalisedLBO(molecule):
    """returns a dict of all atoms with 
    delocalised bonds and their (delocalized-only) VBS"""
    
    def TerminusCounter(atomlist):
        """Counts the number of termini in the delocalized bond system.
        takes a list of all atoms and list of all bonds in the delocalized system as inputs."""
        NTerminus = 0
        for member in atomlist:
            connectivity = 0
            for bond in member.bonds:
                if bond.bond_type == 'Delocalised':
                    connectivity += 1
            if connectivity is 1:
                NTerminus += 1
        return(NTerminus)
    
    def delocal_crawl (atomlist):
        """returns a list of all atoms in a delocalised system,
        takes at least one delocalised bonding atom as input"""
        for delocatom in atomlist:
            for bond in delocatom.bonds:
                if bond.bond_type == 'Delocalised':
                    for member in bond.atoms:
                        if not member in atomlist:
                            atomlist.append(member)
                            return(delocal_crawl(atomlist))
        return(atomlist)
    
    delocal_dict = {}
    for atom in molecule.atoms:
        if all([
            (any(
                bond.bond_type == 'Delocalised'
                for bond in atom.bonds
                )),
                (not atom in delocal_dict)
            ]):
            delocal_dict[atom] = []
            delocal_system = delocal_crawl([atom])
            NTerminus = TerminusCounter(delocal_system)
            for datom in delocal_system:
                connectivity = 0
                delocLBO = 0
                for neighbour in datom.neighbours:
                    if neighbour in delocal_system:
                        connectivity += 1
                if connectivity == 1:
                    #terminus
                    delocLBO = (NTerminus+1)/NTerminus
                if connectivity > 1 :
                    #node
                    delocLBO = (connectivity+1)/connectivity
                delocal_dict[datom] = delocLBO
    return(delocal_dict)

def iVBS_FormalCharge(atom):
    """determines the formal charge on an atom
    that is NOT part of an aromatic or delocalized
    bond system"""
    VBO = 0
    if atom.is_metal:
        return(VBO)
    CN = 0
    for neighbour in atom.neighbours:
        if not neighbour.is_metal:
            CN+=1
    valence = valence_e(atom)
    charge = 0
    for bond in atom.bonds:
        if any(batom.is_metal
        for batom in bond.atoms):
            VBO += 0
        #Each bond contributes to Ligand Bond Order according to its type
        elif bond.bond_type == 'Single':
            VBO +=1
        elif bond.bond_type == 'Double':
            VBO +=2
        elif bond.bond_type == 'Triple':
            VBO +=3
        elif bond.bond_type == 'Quadruple':
            VBO +=4
    #need the unpaired electrons
    unpaired_e = (4 - abs(4 - valence))
    #expanded valences require special handling
    if VBO <= (unpaired_e):
        charge = VBO - unpaired_e
    # Expanded (2e) valences:
    elif (VBO > unpaired_e) and (VBO < valence):
        diff = VBO - unpaired_e
        if diff <= 2:
            UPE = valence - unpaired_e - 2
        elif diff <= 4:
            UPE = valence - unpaired_e - 4
        elif diff <= 6:
            UPE = valence - unpaired_e - 6
        elif diff <= 8:
            UPE = valence - unpaired_e - 8
        charge = valence - (VBO + UPE)
    elif VBO >= (valence):
        charge = valence - VBO
    return(charge)

def get_CN(atom):
    """returns the coordination number of an atom"""
    CN = 0
    for neighbour in atom.neighbours:
        if not neighbour.is_metal:
            CN+=1
    return(CN)

def valence_e(elmnt):
    """returns the number of valence electrons
    of an element"""
    atom = mendeleev.element(elmnt.atomic_symbol)
    if atom.block == 's':
        valence = atom.group_id
    if atom.block == 'p':
        valence = (atom.group_id - 10)
    if atom.block == 'd':
        valence = atom.group_id
    if atom.block == 'f':
        if atom.atomic_number in range(56,72):
            valence = atom.atomic_number - 57 + 3
        elif atom.atomic_number in range(88,104):
            valence = atom.atomic_number - 89 + 3
    if atom.group_id == 18:
        valence = 8
    if atom.symbol == "He":
        valence = 2
    return valence

def carbocation_check (atom):
    """geometry checker for carbocations/carbanions returns
    tetrahedral(anion) or trigonal(cation) depending on bond angles"""
    abc = []
    #get atom neighbours
    for neighbours in atom.neighbours:
        if not neighbours.is_metal:
            abc.append(neighbours)
    #get all three relevant bond angles
    angle1 = descriptors.MolecularDescriptors.atom_angle(abc[0], atom, abc[1])
    angle2 = descriptors.MolecularDescriptors.atom_angle(abc[0], atom, abc[2])
    angle3 = descriptors.MolecularDescriptors.atom_angle(abc[1], atom, abc[2])
    #average the angels
    AVGangle = abs(angle1 + angle2 + angle3)/3
    #take the difference between the averaged bond angles and 
    #ideal trigonal planar/tetrahedral bond angles
    tet = abs(AVGangle - 109.5)
    trig = abs(AVGangle - 120)
    if tet < trig:
        return('tetrahedral')
    if trig < tet:
        return('trigonal')

def carbene_type(atom):
    """distinguishes between singlet and triplet carbenes,
    only input suspected carbene atoms (2-coordinate carbon II)"""
    #get alpha-atoms
    alpha = atom.neighbours
    alpha_type = []
    #get element symbols for alpha atoms
    for a in alpha:
        if not a.is_metal:
            alpha_type.append(a.atomic_symbol)
    # if any alpha atom is a heteroatom, return "singlet"
    # these are Fischer carbenes
    for a in alpha_type:
        if not any([(a == 'C'),
                    (a == 'H')]):
            return('singlet')
    # if the carbene C is in a heterocycle,
    # return "singlet"
    # there are Arduengo carbenes (NHCs, CAACs)
    if atom.is_cyclic == True:
        for ring in atom.rings:
            for species in ring.atoms:
                if not species.atomic_symbol == 'C':
                    return('singlet')
    # for all other carbenes, return "triplet"
    # these are Schrock carbenes
    return('triplet')    

def hapticity(atom, metalsite):
    """returns true if a ligand binding site
    is hapto. Otherwise false"""
    for atom2 in atom.neighbours:
        if not atom2.is_metal:
            if any(
            n2.label == metalsite.label for n2 in atom2.neighbours
            ):
                return(True)
    return(False)

def bridging(atom):
    """Returns the number of metal atoms that 
    a binding site bridges"""
    bridge = 0
    for n in atom.neighbours:
        if n.is_metal:
            bridge += 1
    return (bridge)

def iVBS_Oxidation_Contrib(unique_atoms, rVBO, dVBO):
    """determines the oxidation state contribution for all
    unique atoms in a MOF. Returns a dictionary of Atom:Oxidaiton_Contribution
    pairs. Takes the unique sites in the MOF (without metal), the MOF molecule
    object (without metal) and a list of metal-binding sites as inputs"""
    VBS = 0
    CN = 0
    valence = 0
    oxi_contrib = {}
    # for each unique atom
    for atom in unique_atoms:
        # assign valence-bond-sum
        VBS = assign_VBS(atom, rVBO, dVBO)
        #determine coordination number
        CN = get_CN(atom)
        #  determine number of valence electrons
        valence = valence_e(atom)
        # get number of unpaired electrons in the free element
        unpaired_e = (4 - abs(4 - valence))

        #  metals do not contribute:
        if  atom.is_metal:
            oxi_contrib[atom] = 0
        # Normal valences:
        elif VBS <= (unpaired_e):
            oxi_contrib[atom] = unpaired_e - VBS
        # Expanded (2e) valences:
        elif (VBS > unpaired_e) and (VBS < valence):
            diff = VBS - unpaired_e
            if diff <= 2:
                UPE = valence - unpaired_e - 2
            elif diff <= 4:
                UPE = valence - unpaired_e - 4
            elif diff <= 6:
                UPE = valence - unpaired_e - 6
            elif diff <= 8:
                UPE = valence - unpaired_e - 8
            oxi_contrib[atom] = VBS + UPE - valence
        elif VBS >= (valence):
            oxi_contrib[atom] = VBS - valence
        
        # need to check for 3-coordinate carbocations,
        # 3-coordinate carbanions, carbenes, and heavier
        # homologues (these are not immediately detectable)
        if any([
            (atom.atomic_symbol == 'C'),
            (atom.atomic_symbol == 'Si'),
            (atom.atomic_symbol == 'Ge'),
            (atom.atomic_symbol is 'Pb')]):
             if not atom in rVBO:
              # 3 coordinate and VBS 3 could be
              # carbanion or carbocation
              if VBS == 3 and CN == 3:
                  geom = carbocation_check(atom)
                  if geom == "trigonal":
                      oxi_contrib[atom] = -1
                  if geom == "tetrahedral":
                      oxi_contrib[atom] = 1
             # VBS 2 and 2 coordinate is carbene,
             # but singlet or triplet?
             if VBS == 2 and CN == 2:
               carbene = carbene_type(atom)
               if carbene == 'singlet':
                 oxi_contrib[atom] = 2
               if carbene == 'triplet':
                 oxi_contrib[atom] = 0

        # Nitro groups frequently have both N-O bonds assigned
        # as double bonds, giving incorrect VBS of 5
        # and oxidation contribution of -2
        #this block catches this and applies a fix
        if all([
            (atom.atomic_symbol == 'N'),
            (VBS == 5 and CN == 3),
        ]):
            N_sphere1 = atom.neighbours
            O_count = 0
            for neighbour in N_sphere1:
                if neighbour.atomic_symbol == 'O':
                    O_count += 1
            geom = carbocation_check(atom)
            if O_count == 2 and geom == "trigonal":
                oxi_contrib[atom] = 0
               
    return (oxi_contrib)

def redundantAON(AON, molecule):
    redAON = {}
    for rsite1 in molecule.atoms:
        for usite1 in AON:
            redAON[usite1] = AON[usite1]
            if rsite1.label == usite1.label:
                redAON[rsite1] = AON[usite1]
    return(redAON)

def binding_domain(binding_sites, AON, molecule, usites):
    """Due to geometric distortions (as expected in a MOF),
    delocalized bonds are not always correctly assigned
    This function will build 'bonding domains' to fix this 
    (see methodology section of paper)"""
    def arom_domains(site, usites, aromlist, bondset):
        """recursive crawler makes aromatic domains"""
        for bond in site.bonds:
            bondset.add(bond)
        for bond in bondset:
            for member in bond.atoms:
                if(
                    all([
                        (not member in aromlist),
                        (not member.is_metal),
                        (any(
                            mbond.bond_type == 'Aromatic'
                            for mbond in member.bonds))
                        ])
                    ):
                    aromlist.append(member)
                    for mbond in member.bonds:
                        bondset.add(mbond)
                    return(arom_domains(site, usites, aromlist, bondset))
        #aromlist currently contains non-unique instances of atoms
        #this will cause problems further down the line, so correct
        for index, member in enumerate(aromlist):
            aromlist[index] = usites[member.label]
        return(aromlist)
    def deloc_domains(site, usites, AON, molecule, deloclist, bondset, checked_bonds):
        """recursive crawler makes domains"""
        for bond in site.bonds:
            if not bond in bondset:
                bondset.add(bond)
        for bond in bondset:
            if not bond in checked_bonds:
                for member in bond.atoms:
                    if( 
                    all([
                    (not member in deloclist),
                    (not member.is_metal),
                    (not any(mbond.bond_type == 'Aromatic'
                    for mbond in member.bonds)),
                    (any([
                        (len(molecule.shortest_path_bonds(site, member)) <= 2),
                        (bond.bond_type == 'Delocalised'),
                        (bond.is_conjugated),
                        (all([
                            (bond.bond_type == 'Single'),
                            (not AON[member] == 0),
                            ])),
                        (all([
                            (not any(mbond.bond_type == 'Single'
                            for mbond in member.bonds)),
                            (not any(mbond.bond_type == 'Aromatic'
                            for mbond in member.bonds)),
                            (not any(mbond.bond_type == 'Delocalised'
                            for mbond in member.bonds))
                            ]))
                        ]))
                        ])
                        ):
                        deloclist.append(member)
                        for mbond in member.bonds:
                            bondset.add(mbond)
                checked_bonds.add(bond)
                return(deloc_domains(site, usites, AON, molecule, deloclist, bondset, checked_bonds))
        #deloclist currently contains non-unique instances of atoms
        #this will cause problems further down the line, so correct
        for index, member in enumerate(deloclist):
            deloclist[index] = usites[member.label]
        return(deloclist)
    sitedomain={}
    for site in binding_sites:
        if not site.is_metal == True:
            if any(sbond.bond_type == 'Aromatic' for sbond in site.bonds):
                sitedomain[site] = arom_domains(site, usites, aromlist=[site], bondset=set())
            if not any(sbond.bond_type == 'Aromatic' for sbond in site.bonds):
                sitedomain[site] = deloc_domains(site, usites, AON, molecule, deloclist=[site], bondset=set(), checked_bonds=set())
    
    for site in sitedomain:
        olapset = set()
        for site2 in sitedomain:
            for member in sitedomain[site]:
                if member in sitedomain[site2]:
                    olapset.add(site2)
        for olap in olapset:
            sitedomain[site] = list(set(sitedomain[site]) | set(sitedomain[olap]))
            sitedomain[olap] = sitedomain[site]
    return(sitedomain)

def binding_contrib(binding_sphere, binding_sites, AON):
    """Within a binding domain, oxidation state contribution
    is assumed to be equally distributed across all binding sites"""
    site_contrib = {}
    for site in binding_sphere:
        site_contrib[site] = 0
        nbinding = 0
        for member in binding_sphere[site]:
            if member in binding_sites:
                nbinding += 1
            site_contrib[site] += AON[member]
        site_contrib[site] /= nbinding
    return (site_contrib)

def outer_sphere_domain (uniquesites,binding_domains):
    outer_sphere = []
    for site in uniquesites:
        if (all([
            (not any(site in binding_domains[domain] 
                for domain in binding_domains)),
            (not site.is_metal),
            ])):
            outer_sphere.append(site)
    return(outer_sphere)

def outer_sphere_contrib (outer_sphere, AON):
    contrib = 0
    for site in outer_sphere:
        contrib += AON[site]
    return(contrib)

def get_metal_networks(ligand_sites, binding_sphere, bindingAON):
    def network_crawl(ligand_sites, binding_sphere, bindingAON, metal_networks, checked_sites, group):
        for metal in group:
            #This block will find all metals connected to an input metal by metal-metal bonds
            checked_sites.append(metal)
            for neighbour in metal.neighbours:
                if neighbour.is_metal:
                    if not neighbour in checked_sites:
                        checked_sites.append(neighbour)
                        for site in ligand_sites:
                            if neighbour.label == site.label:
                                if not site in group:
                                    group.append(site)
                        return network_crawl(ligand_sites, binding_sphere, bindingAON, metal_networks, checked_sites, group)
            # this block will find all metals connected to an input metal by
            # conjugation and delocalized charge ligands
            # metals connected through NEUTRAL ligands will be ignored
            for site in ligand_sites[metal]:
                if all([
                    (not bindingAON[site] == 0),
                    (not site in checked_sites)
                ]):
                    for dsite in binding_sphere[site]:
                        if all([
                            (not dsite in checked_sites),
                            (dsite in binding_sphere)
                        ]):
                            checked_sites.append(dsite)
                            for environ in dsite.neighbours:
                                if environ.is_metal:
                                    if environ in ligand_sites:
                                        if all([                                            
                                            (all(not environ in network for network in metal_networks)),
                                            (not environ in group)
                                        ]):
                                            group.append(environ)
                                    else:
                                        for umetal in ligand_sites:
                                            if all([
                                                (umetal.label == environ.label),
                                                (all(not umetal in network for network in metal_networks)),
                                                (not umetal in group)
                                            ]):
                                                group.append(umetal)
                    return network_crawl(ligand_sites, binding_sphere, bindingAON, metal_networks, checked_sites, group)
        return (group)

    metal_networks = []
    for metal in ligand_sites:
        if (all(not metal in network for network in metal_networks)):
            metal_networks.append(
                network_crawl(ligand_sites,
                                          binding_sphere,
                                          bindingAON, 
                                          metal_networks, 
                                          checked_sites=[], 
                                          group=[metal]))

    network_dict = {}
    for network in metal_networks:
        for metal in network:
            network_dict[metal] = network
    return (network_dict)

def distribute_ONEC(sONEC, metal_networks, IEs, ONP, highest_known_ON, metal_CN):
    """function to redistribute electron count and oxidation state
    returns an updated dictionary in the format metal:[ON,EC]"""
    
    def recursive_distributor_single_network(iONEC, available_charge, sorted_metals, IEs, ONP, highest_known_ON):
        """run this function after tallying available network charge and sorting network metals by
        element type. Will distribute network charge according to ionization energy and probability until all
        charge is distributed."""
        #initialize working dictionary
        dONEC = {}
        dONEC = dict(iONEC)
        
        # positive contribution?
        if available_charge > 0:
            
            #get list of improbable and improbable next oxidations
            prob_metal_type = []
            improb_metal_type = []
            for metal_type in sorted_metals:
                try:
                    prob = float(100*ONP[metal_type][math.floor(dONEC[sorted_metals[metal_type][0]][0])+1])
                except IndexError:
                    prob = 0
                if prob >= 10:
                    prob_metal_type.append(metal_type)
                else:
                    improb_metal_type.append(metal_type)
            
            #if only one metal type has a probable next oxidation state, do that
            if len(prob_metal_type) == 1:
                lowestMetal = prob_metal_type[0]

            #if more than one metal type has a probable next oxidation state,
            #determine next lowest ionization energy among probable next
            #oxidation states
            elif len(prob_metal_type) > 1:      
                #find lowest next ionization energy
                for metal_type in prob_metal_type:
                    #metal in a negative oxidation state? Use IE = 0.
                    if dONEC[sorted_metals[metal_type][0]][0] < 0:
                        currentIE = 0
                    #metal oxidation state at or higher than highest known? Set IE arbitrarily high.
                    elif dONEC[sorted_metals[metal_type][0]][0] >= highest_known_ON[metal_type]:
                        currentIE = 9999
                    #otherwise, use the appropriate IE.
                    else:
                        currentIE = float(IEs[metal_type][math.floor(dONEC[sorted_metals[metal_type][0]][0])])
                    if not 'lowestIE' in locals():
                        lowestIE = currentIE
                        lowestMetal = metal_type
                    else:
                        if currentIE < lowestIE:
                            lowestIE = currentIE
                            lowestMetal = metal_type 

            #if there is no probable next oxidation state available,
            #determine lowest ionization energy among improbable next oxidation states
            elif len(prob_metal_type) == 0:
                #find lowest next ionization energy
                for metal_type in improb_metal_type:
                    #metal in a negative oxidation state? Use IE = 0.
                    if dONEC[sorted_metals[metal_type][0]][0] < 0:
                        currentIE = 0
                    #metal oxidation state at or higher than highest known? Set IE arbitrarily high.
                    elif dONEC[sorted_metals[metal_type][0]][0] >= highest_known_ON[metal_type]:
                        currentIE = 9999
                    #otherwise, use the appropriate IE.
                    else:
                        currentIE = float(IEs[metal_type][math.floor(dONEC[sorted_metals[metal_type][0]][0])])
                    if not 'lowestIE' in locals():
                        lowestIE = currentIE
                        lowestMetal = metal_type
                    else:
                        if currentIE < lowestIE:
                            lowestIE = currentIE
                            lowestMetal = metal_type

            #distribute one ionization energy level worth of charge
            if available_charge >= len(sorted_metals[lowestMetal]):
                for metal in sorted_metals[lowestMetal]:
                    dONEC[metal][0] += 1
                    available_charge -= 1
            elif available_charge < len(sorted_metals[lowestMetal]):
                for metal in sorted_metals[lowestMetal]:
                    dONEC[metal][0] += available_charge/(len(sorted_metals[lowestMetal]))
                available_charge = 0

        # negative contribution?
        if available_charge < 0:
            #get list of improbable and improbable next oxidations
            prob_metal_type = []
            improb_metal_type = []
            for metal_type in sorted_metals:
                try:
                    prob = float(100*ONP[metal_type][math.floor(dONEC[sorted_metals[metal_type][0]][0])-1])
                except IndexError:
                    prob = 0
                if prob >= 10:
                    prob_metal_type.append(metal_type)
                else:
                    improb_metal_type.append(metal_type)

            #if only one metal type has a probable next oxidation state, do that
            if len(prob_metal_type) == 1:
                highestMetal = prob_metal_type[0]
            
            #if more than one metal type has a probable next oxidation state,
            #determine next highest ionization energy among probable next
            #oxidation states
            elif len(prob_metal_type) > 1:
                for metal_type in prob_metal_type:
                    #metal in a negative oxidation state? Use IE = 0.
                    if dONEC[sorted_metals[metal_type][0]][0] <= 0:
                        currentIE = 0
                    #metal oxidation state higher than highest known? Set IE arbitrarily high.
                    elif dONEC[sorted_metals[metal_type][0]][0] > highest_known_ON[metal_type]:
                        currentIE = 9999
                    #otherwise, use the appropriate IE.
                    else:
                        if not any([
                        (math.isclosey((dONEC[sorted_metals[metal_type][0]][0] % 1), 0, abs_tol=0.0001)),
                        (math.isclose((dONEC[sorted_metals[metal_type][0]][0] % 1), 1, abs_tol=0.0001))
                        ]):
                            currentIE = float(IEs[metal_type][math.floor(dONEC[sorted_metals[metal_type][0]][0])])
                        else:
                            currentIE = float(IEs[metal_type][math.floor(dONEC[sorted_metals[metal_type][0]][0])-1])
                    if not 'highestIE' in locals():
                        highestIE = currentIE
                        highestMetal = metal_type
                    else:
                        if currentIE > highestIE:
                            highestIE = currentIE
                            highestMetal = metal_type    

            #if no probable next oxidation states are available,
            #determine next highest ionization energy among probable next
            #oxidation states
            elif len(improb_metal_type) > 0:
                for metal_type in improb_metal_type:
                    #metal in a negative oxidation state? Use IE = 0.
                    if dONEC[sorted_metals[metal_type][0]][0] <= 0:
                        currentIE = 0
                    #metal oxidation state higher than highest known? Set IE arbitrarily high.
                    elif dONEC[sorted_metals[metal_type][0]][0] > highest_known_ON[metal_type]:
                        currentIE = 9999
                    #otherwise, use the appropriate IE.
                    else:
                        if not any([
                        (math.isclose((dONEC[sorted_metals[metal_type][0]][0] % 1), 0, abs_tol=0.0001)),
                        (math.isclose((dONEC[sorted_metals[metal_type][0]][0] % 1), 1, abs_tol=0.0001))
                        ]):
                            currentIE = float(IEs[metal_type][math.floor(dONEC[sorted_metals[metal_type][0]][0])])
                        else:
                            currentIE = float(IEs[metal_type][math.floor(dONEC[sorted_metals[metal_type][0]][0])-1])
                    if not 'highestIE' in locals():
                        highestIE = currentIE
                        highestMetal = metal_type
                    else:
                        if currentIE > highestIE:
                            highestIE = currentIE
                            highestMetal = metal_type

            #distribute one ionization energy level worth of charge
            if (-1*available_charge) >= len(sorted_metals[highestMetal]):
                for metal in sorted_metals[highestMetal]:
                    dONEC[metal][0] -= 1
                    available_charge += 1
            elif (-1*available_charge) < len(sorted_metals[highestMetal]):
                for metal in sorted_metals[highestMetal]:
                    dONEC[metal][0] += available_charge/(len(sorted_metals[highestMetal]))
                available_charge = 0
        
        #if all charge has been distributed, we're done, otherwise, roll again
        if available_charge == 0:
            return(dONEC)
        else:
            return (recursive_distributor_single_network(dONEC, available_charge, sorted_metals, IEs, ONP, highest_known_ON))
    
    #operate on each network individually
    distributed_ONEC = {}
    for network in metal_networks: 

        #sort metals by element type
        sorted_metals = {}
        for metal in metal_networks[network]:
            sorted_metals[metal.atomic_symbol] = []
        for metal in metal_networks[network]:
            sorted_metals[metal.atomic_symbol].append(metal)
       
        # tally up network charge to be distributed
        # and initialize distributed charge for network metals
        network_charge = 0
        for metal in metal_networks[network]:
            network_charge += sONEC[metal][0]
            distributed_ONEC[metal] = [0]
            distributed_ONEC[metal].append(int(sONEC[metal][1]))
        
        #recursively distribute network charge according to ionization energy
        distributed_ONEC = recursive_distributor_single_network(distributed_ONEC, network_charge, sorted_metals, IEs, ONP, highest_known_ON)
        
        # finally, adjust electron count to new oxidation state (OiL RiG)
        for metal in metal_networks[network]:
            distributed_ONEC[metal][1] = (valence_e(metal) + (2 * metal_CN[metal]) - distributed_ONEC[metal][0])

    return (distributed_ONEC)

def distribute_OuterSphere(sONEC, outer_sphere_charge, IEs, ONP, highest_known_ON, metal_CN):
    """function to redistribute electron count and oxidation state
    returns an updated dictionary in the format metal:[ON,EC]"""
    
    def recursive_distributor(iONEC, available_charge, IEs, ONP, highest_known_ON):
        """run this function after tallying available network charge and sorting network metals by
        element type. Will distribute network charge according to ionization energy and highest allowable oxidation state until all
        charge is distributed."""
        #initialize working dictionary
        dONEC = {}
        dONEC = dict(iONEC)
        
        # positive contribution?
        if available_charge > 0:
            #get list of improbable and improbable next oxidations
            prob_metals = []
            improb_metals = []
            for metal in dONEC:
                try:
                    prob = float(100*ONP[metal.atomic_symbol][math.floor(dONEC[metal][0])+1])
                except IndexError:
                    prob = 0
                if prob >= 10:
                    prob_metals.append(metal)
                else:
                    improb_metals.append(metal)
            
            if len(prob_metals) == 1:
                lowestMetals = prob_metals
            elif len(prob_metals) > 1:
                for metal in prob_metals:
                    #metal in a negative oxidation state? Use IE = 0.
                    if dONEC[metal][0] < 0:
                        currentIE = 0
                    #metal oxidation state at or higher than highest known? Set IE arbitrarily high.
                    elif dONEC[metal][0] >= highest_known_ON[metal.atomic_symbol]:
                        currentIE = 9999
                    #otherwise, use the appropriate IE.
                    else:
                        currentIE = float(IEs[metal.atomic_symbol][math.floor(dONEC[metal][0])])
                    if not 'lowestIE' in locals():
                        lowestIE = currentIE
                        lowestMetals = [metal]
                    else:
                        if currentIE == lowestIE:
                            lowestMetals.append(metal)
                        if currentIE < lowestIE:
                            lowestIE = currentIE
                            lowestMetals = [metal]  
            elif len(prob_metals) == 0:
                for metal in improb_metals:
                    #metal in a negative oxidation state? Use IE = 0.
                    if dONEC[metal][0] < 0:
                        currentIE = 0
                    #metal oxidation state at or higher than highest known? Set IE arbitrarily high.
                    elif dONEC[metal][0] >= highest_known_ON[metal.atomic_symbol]:
                        currentIE = 9999
                    #otherwise, use the appropriate IE.
                    else:
                        currentIE = float(IEs[metal.atomic_symbol][math.floor(dONEC[metal][0])])
                    if not 'lowestIE' in locals():
                        lowestIE = currentIE
                        lowestMetals = [metal]
                    else:
                        if currentIE == lowestIE:
                            lowestMetals.append(metal)
                        if currentIE < lowestIE:
                            lowestIE = currentIE
                            lowestMetals = [metal] 
            #distribute one ionization energy level worth of charge
            if available_charge >= len(lowestMetals):
                for metal in lowestMetals:
                    dONEC[metal][0] += 1
                    available_charge -= 1
            elif available_charge < len(lowestMetals):
                for metal in lowestMetals:
                    dONEC[metal][0] += available_charge/(len(lowestMetals))
                available_charge = 0

        # negative contribution?
        if available_charge < 0:
            #get list of improbable and improbable next oxidations
            prob_metals = []
            improb_metals = []
            for metal in dONEC:
                try:
                    prob = float(100*ONP[metal.atomic_symbol][math.floor(dONEC[metal][0])-1])
                except IndexError:
                    prob = 0
                if prob >= 10:
                    prob_metals.append(metal)
                else:
                    improb_metals.append(metal)
            
            if len(prob_metals) == 1:
                highestMetals = prob_metals
            elif len(prob_metals) > 1:
                for metal in prob_metals:
                    #metal in a negative oxidation state? Use IE = 0.
                    if dONEC[metal][0] <= 0:
                        currentIE = 0
                    #metal oxidation state higher than highest known? Set IE arbitrarily high.
                    elif dONEC[metal][0] > highest_known_ON[metal.atomic_symbol]:
                        currentIE = 9999
                    #otherwise, use the appropriate IE.
                    else:
                        if not any([
                        (math.isclose((dONEC[metal][0] % 1), 0, abs_tol=0.0001)),
                        (math.isclose((dONEC[metal][0] % 1), 1, abs_tol=0.0001))
                        ]):
                            currentIE = float(IEs[metal.atomic_symbol][math.floor(dONEC[metal][0])])
                        else:
                            currentIE = float(IEs[metal.atomic_symbol][math.floor(dONEC[metal][0])-1])
                    if not 'highestIE' in locals():
                        highestIE = currentIE
                        highestMetals = [metal]
                    else:
                        if currentIE == highestIE:
                            highestMetals.append(metal)
                        if currentIE > highestIE:
                            highestIE = currentIE
                            highestMetals = [metal]
            
            else:
                for metal in improb_metals:
                    #metal in a negative oxidation state? Use IE = 0.
                    if dONEC[metal][0] <= 0:
                        currentIE = 0
                    #metal oxidation state higher than highest known? Set IE arbitrarily high.
                    elif dONEC[metal][0] > highest_known_ON[metal.atomic_symbol]:
                        currentIE = 9999
                    #otherwise, use the appropriate IE.
                    else:
                        if not any([
                        (math.isclose((dONEC[metal][0] % 1), 0, abs_tol=0.0001)),
                        (math.isclose((dONEC[metal][0] % 1), 1, abs_tol=0.0001))
                        ]):
                            currentIE = float(IEs[metal.atomic_symbol][math.floor(dONEC[metal][0])])
                        else:
                            currentIE = float(IEs[metal.atomic_symbol][math.floor(dONEC[metal][0])-1])
                    if not 'highestIE' in locals():
                        highestIE = currentIE
                        highestMetals = [metal]
                    else:
                        if currentIE == highestIE:
                            highestMetals.append(metal)
                        if currentIE > highestIE:
                            highestIE = currentIE
                            highestMetals = [metal]

            #distribute one ionization energy level worth of charge
            if (-1*available_charge) >= len(highestMetals):
                for metal in highestMetals:
                    dONEC[metal][0] -= 1
                    available_charge += 1
            elif (-1*available_charge) < len(highestMetals):
                for metal in highestMetals:
                    dONEC[metal][0] += available_charge/(len(highestMetals))
                available_charge = 0
        
        #if all charge has been distributed, we're done, otherwise, roll again
        if available_charge == 0:
            return(dONEC)
        else:
            return (recursive_distributor(dONEC, available_charge, IEs, ONP, highest_known_ON))

    if outer_sphere_charge == 0:
        return (sONEC)
    
    distributed_ONEC = {}
    
    # initialize dictionary for charge distribution
    for metal in sONEC:
        distributed_ONEC[metal] = [sONEC[metal][0]]
        distributed_ONEC[metal].append(int(sONEC[metal][1]))
        
    #recursively distribute network charge according to ionization energy
    distributed_ONEC = recursive_distributor(distributed_ONEC, outer_sphere_charge, IEs, ONP, highest_known_ON)
        
    # finally, adjust electron count to new oxidation state (OiL RiG)
    for metal in sONEC:
        distributed_ONEC[metal][1] = (valence_e(metal) + (2 * metal_CN[metal]) - distributed_ONEC[metal][0])

    return (distributed_ONEC)

def global_charge_distribution(metalONdict, IEs, ONP, highest_known_ON, metal_CN):
    global_ONEC = {}
    
    def recursive_distributor_global(iONEC, available_charge, sorted_metals, IEs, ONP, highest_known_ON):
        """run this function after tallying available network charge and sorting network metals by
        element type. Will distribute network charge according to ionization energy and probability until all
        charge is distributed."""
        #initialize working dictionary
        dONEC = {}
        dONEC = dict(iONEC)
        
        # positive contribution?
        if available_charge > 0:
            
            #get list of improbable and improbable next oxidations
            prob_metal_type = []
            improb_metal_type = []
            for metal_type in sorted_metals:
                try:
                    prob = float(100*ONP[metal_type][math.floor(dONEC[sorted_metals[metal_type][0]][0])+1])
                except IndexError:
                    prob = 0
                if prob >= 10:
                    prob_metal_type.append(metal_type)
                else:
                    improb_metal_type.append(metal_type)
            
            #if only one metal type has a probable next oxidation state, do that
            if len(prob_metal_type) == 1:
                lowestMetal = prob_metal_type[0]

            #if more than one metal type has a probable next oxidation state,
            #determine next lowest ionization energy among probable next
            #oxidation states
            elif len(prob_metal_type) > 1:      
                #find lowest next ionization energy
                for metal_type in prob_metal_type:
                    #metal in a negative oxidation state? Use IE = 0.
                    if dONEC[sorted_metals[metal_type][0]][0] < 0:
                        currentIE = 0
                    #metal oxidation state at or higher than highest known? Set IE arbitrarily high.
                    elif dONEC[sorted_metals[metal_type][0]][0] >= highest_known_ON[metal_type]:
                        currentIE = 9999
                    #otherwise, use the appropriate IE.
                    else:
                        currentIE = float(IEs[metal_type][math.floor(dONEC[sorted_metals[metal_type][0]][0])])
                    if not 'lowestIE' in locals():
                        lowestIE = currentIE
                        lowestMetal = metal_type
                    else:
                        if currentIE < lowestIE:
                            lowestIE = currentIE
                            lowestMetal = metal_type 

            #if there is no probable next oxidation state available,
            #determine lowest ionization energy among improbable next oxidation states
            elif len(prob_metal_type) == 0:
                #find lowest next ionization energy
                for metal_type in improb_metal_type:
                    #metal in a negative oxidation state? Use IE = 0.
                    if dONEC[sorted_metals[metal_type][0]][0] < 0:
                        currentIE = 0
                    #metal oxidation state at or higher than highest known? Set IE arbitrarily high.
                    elif dONEC[sorted_metals[metal_type][0]][0] >= highest_known_ON[metal_type]:
                        currentIE = 9999
                    #otherwise, use the appropriate IE.
                    else:
                        currentIE = float(IEs[metal_type][math.floor(dONEC[sorted_metals[metal_type][0]][0])])
                    if not 'lowestIE' in locals():
                        lowestIE = currentIE
                        lowestMetal = metal_type
                    else:
                        if currentIE < lowestIE:
                            lowestIE = currentIE
                            lowestMetal = metal_type

            #distribute one ionization energy level worth of charge
            if available_charge >= len(sorted_metals[lowestMetal]):
                for metal in sorted_metals[lowestMetal]:
                    dONEC[metal][0] += 1
                    available_charge -= 1
            elif available_charge < len(sorted_metals[lowestMetal]):
                for metal in sorted_metals[lowestMetal]:
                    dONEC[metal][0] += available_charge/(len(sorted_metals[lowestMetal]))
                available_charge = 0

        # negative contribution?
        if available_charge < 0:
            #get list of improbable and improbable next oxidations
            prob_metal_type = []
            improb_metal_type = []
            for metal_type in sorted_metals:
                try:
                    prob = float(ONP[metal_type][math.floor(dONEC[sorted_metals[metal_type][0]][0])-1])
                except IndexError:
                    prob = 0
                if prob >= 10:
                    prob_metal_type.append(metal_type)
                else:
                    improb_metal_type.append(metal_type)

            #if only one metal type has a probable next oxidation state, do that
            if len(prob_metal_type) == 1:
                highestMetal = prob_metal_type[0]
            
            #if more than one metal type has a probable next oxidation state,
            #determine next highest ionization energy among probable next
            #oxidation states
            elif len(prob_metal_type) > 1:
                for metal_type in prob_metal_type:
                    #metal in a negative oxidation state? Use IE = 0.
                    if dONEC[sorted_metals[metal_type][0]][0] <= 0:
                        currentIE = 0
                    #metal oxidation state higher than highest known? Set IE arbitrarily high.
                    elif dONEC[sorted_metals[metal_type][0]][0] > highest_known_ON[metal_type]:
                        currentIE = 9999
                    #otherwise, use the appropriate IE.
                    else:
                        if not any([
                        (math.isclosey((dONEC[sorted_metals[metal_type][0]][0] % 1), 0, abs_tol=0.0001)),
                        (math.isclose((dONEC[sorted_metals[metal_type][0]][0] % 1), 1, abs_tol=0.0001))
                        ]):
                            currentIE = float(IEs[metal_type][math.floor(dONEC[sorted_metals[metal_type][0]][0])])
                        else:
                            currentIE = float(IEs[metal_type][math.floor(dONEC[sorted_metals[metal_type][0]][0])-1])
                    if not 'highestIE' in locals():
                        highestIE = currentIE
                        highestMetal = metal_type
                    else:
                        if currentIE > highestIE:
                            highestIE = currentIE
                            highestMetal = metal_type    

            #if no probable next oxidation states are available,
            #determine next highest ionization energy among probable next
            #oxidation states
            elif len(improb_metal_type) > 0:
                for metal_type in improb_metal_type:
                    #metal in a negative oxidation state? Use IE = 0.
                    if dONEC[sorted_metals[metal_type][0]][0] <= 0:
                        currentIE = 0
                    #metal oxidation state higher than highest known? Set IE arbitrarily high.
                    elif dONEC[sorted_metals[metal_type][0]][0] > highest_known_ON[metal_type]:
                        currentIE = 9999
                    #otherwise, use the appropriate IE.
                    else:
                        if not any([
                        (math.isclose((dONEC[sorted_metals[metal_type][0]][0] % 1), 0, abs_tol=0.0001)),
                        (math.isclose((dONEC[sorted_metals[metal_type][0]][0] % 1), 1, abs_tol=0.0001))
                        ]):
                            currentIE = float(IEs[metal_type][math.floor(dONEC[sorted_metals[metal_type][0]][0])])
                        else:
                            currentIE = float(IEs[metal_type][math.floor(dONEC[sorted_metals[metal_type][0]][0])-1])
                    if not 'highestIE' in locals():
                        highestIE = currentIE
                        highestMetal = metal_type
                    else:
                        if currentIE > highestIE:
                            highestIE = currentIE
                            highestMetal = metal_type

            #distribute one ionization energy level worth of charge
            if (-1*available_charge) >= len(sorted_metals[highestMetal]):
                for metal in sorted_metals[highestMetal]:
                    dONEC[metal][0] -= 1
                    available_charge += 1
            elif (-1*available_charge) < len(sorted_metals[highestMetal]):
                for metal in sorted_metals[highestMetal]:
                    dONEC[metal][0] += available_charge/(len(sorted_metals[highestMetal]))
                available_charge = 0
        
        #if all charge has been distributed, we're done, otherwise, roll again
        if available_charge == 0:
            return(dONEC)
        else:
            return (recursive_distributor_global(dONEC, available_charge, sorted_metals, IEs, ONP, highest_known_ON))
    
    #sort metals by element type
    sorted_metals = {}
    for metal in metalONdict:
        sorted_metals[metal.atomic_symbol] = []
    for metal in metalONdict:
        sorted_metals[metal.atomic_symbol].append(metal)
       
    # tally up the global charge to be distributed
    # and initialize global ON to 0 for all metals
    global_charge = 0
    for metal in metalONdict:
        global_charge += metalONdict[metal][0]
        global_ONEC[metal] = [0]
        global_ONEC[metal].append(int(metalONdict[metal][1]))
        
    #recursively distribute network charge according to ionization energy
    distributed_ONEC = recursive_distributor_global(global_ONEC, global_charge, sorted_metals, IEs, ONP, highest_known_ON)
        
    # finally, adjust electron count to new oxidation state (OiL RiG)
    for metal in metalONdict:
        global_ONEC[metal][1] = (valence_e(metal) + (2 * metal_CN[metal]) - distributed_ONEC[metal][0])
    return (global_ONEC)
    
def equiv_metals(metalsites, AON, ONEC):
    environ = {}
    for site in metalsites:
        environ[site] = (list())
        for ligand in site.neighbours:
            environ[site].append(f'{ligand.atomic_symbol}{AON[ligand]}')
    equivalent_metals = {}
    metal_groups = {}
    for metal in environ:
        equivalent_metals[metal] = [metal]
        for metal2 in environ:
            if not metal == metal2:
                if all([(metal.atomic_symbol == metal2.atomic_symbol),
                    (math.isclose(ONEC[metal][0],ONEC[metal2][0])),
                    (math.isclose(ONEC[metal][1],ONEC[metal2][1])),
                    (all(i in environ[metal2]
                    for i in environ[metal]))]):
                    equivalent_metals[metal].append(metal2)
    def remove_duplicates(equivalent_metals):
        for metal in equivalent_metals:
            for metal2 in equivalent_metals:
                if not metal == metal2:
                    if all(equivalent_metals[metal][i] in equivalent_metals[metal2] 
                    for i in range(0,len(equivalent_metals[metal]))):
                        equivalent_metals.pop(metal2)
                        return(remove_duplicates(equivalent_metals))
        return(equivalent_metals)
    equivalent_metals = remove_duplicates(equivalent_metals)
    for metal in equivalent_metals:
            metal_groups[ONEC[metal][0]] = equivalent_metals[metal]
    return (metal_groups)

def KnownONs():
    KONs = {}
    with open(os.path.join(cwd,'KnownON.csv')) as ONs:
        for ON in ONs.readlines():
            ONlist = []
            splitON = ON.split(',')
            for split in splitON:
                split.replace(',','')
            while('' in splitON):
                splitON.remove('')
            while('\n' in splitON):
                splitON.remove('\n')
            for i in range(1,len(splitON)):
                ONlist.append(splitON[i])
            KONs[splitON[0]] = ONlist
    ONs.close()
    return(KONs)

def IonizationEnergies():
    KIEs = {}
    IElist = []
    with open(os.path.join(cwd,'Ionization_Energies.csv')) as IEs:
        for IE in IEs.readlines():
            splitIE = IE.split(',')
            for split in splitIE:
                split.replace(',','')
                if(r'\n' in split):
                    split.replace(r'\n', '')
            while('' in splitIE):
                splitIE.remove('')
            while(r'\n' in splitIE):
                splitIE.remove(r'\n')
            IElist.append(splitIE)
        for entry in IElist:
            if entry[1] in KIEs:
                KIEs[(entry[1])].append(entry[2])
            else:
                KIEs[entry[1]] = []
                KIEs[entry[1]].append(entry[2])
    IEs.close()
    return(KIEs)

def HighestKnownONs():
    HKONs = {}
    with open(os.path.join(cwd,'KnownON.csv')) as ONs:
        for ON in ONs.readlines():
            highest = 0
            splitON = ON.split(',')
            for split in splitON:
                split.replace(',','')
            while('' in splitON):
                splitON.remove('')
            while('\n' in splitON):
                splitON.remove('\n')
            for i in range(1,len(splitON)):
                if int(splitON[i]) >= highest:
                    highest = int(splitON[i])
            HKONs[splitON[0]] = int(highest)
    ONs.close()
    return(HKONs)

def ONprobabilities():
    ONP = {}
    with open(os.path.join(cwd,'Oxidation_Probabilities.csv')) as ONPs:
        for ON in ONPs.readlines():
            ONPlist = []
            splitONP = ON.split(',')
            for split in splitONP:
                split.replace(',','')
            while('' in splitONP):
                splitONP.remove('')
            while('\n' in splitONP):
                splitONP.remove('\n')
            for i in range(1,len(splitONP)):
                ONPlist.append(float(splitONP[i]))
            ONP[splitONP[0]] = ONPlist
    ONPs.close()
    return(ONP)
    
def getCN(lsites):
    CNdict = {}
    for metal in lsites:
        CNdict[metal] = 0
        for ligand in lsites[metal]:
            if hapticity(ligand, metal):
                CNdict[metal] += 0.5
            else:
                CNdict[metal] += 1
        for neighbour in metal.neighbours:
            if neighbour.is_metal:
                CNdict[metal] += 0.5
    return(CNdict)

cwd = os.getcwd()
KnownON = KnownONs()
KnownIE = IonizationEnergies()
HighestKnownON = HighestKnownONs()
ONProb = ONprobabilities()
alldata = pd.DataFrame()
ZeroVectorlist = []

for file in os.listdir(cwd):
    if all([
          (any([
              (str(file).endswith('.cif') == True),
              (str(file).endswith('.mol2') == True),
              ])),
          (not str(file).startswith('._')),
          ]):
        print (f'Analyzing {file}...')
        input_cif = (file)
        # read in the .cif, extract the underlying molecule,
        # identify the unique sites, metal sites, binding sites,
        # ect.
        #print ('reading .cif data')
        if str(file).endswith('.cif'):
          try:
            cif = readentry(input_cif)
            mol = cif.molecule
            asymmol = cif.asymmetric_unit_molecule
          except RuntimeError:
            ZeroVectorlist.append(file)
            continue
        if str(file).endswith('.mol2'):
          mol = readSBU(file)
          asymmol = readSBU(file)
        #print ('identifying unique sites')
        uniquesites = get_unique_sites(mol, asymmol)
        usitedict = {}
        for usite in uniquesites:
            usitedict[usite.label] = usite
        #print ('getting metal sites')
        metalsites = get_metal_sites(uniquesites)
        if len(metalsites) == 0:
          nometals = "NO_METALS"
          entrydata = pd.DataFrame({"CIF":[file],
                                         "Metal":[nometals], 
                                         "ON_coordination_ONLY":[nometals],
                                         "EC_coordination_ONLY":[nometals],
                                         "ON_coordination+Outer_Sphere":[nometals],
                                         "EC_coordination+Outer_Sphere":[nometals],
                                         "ON_network_redistributed":[nometals],
                                         "EC_network_redistributed":[nometals],
                                         "ON_network+Outer_Sphere":[nometals],
                                         "EC_network+Outer_Sphere":[nometals],
                                         "Impossible":[nometals],
                                         "Unknown":[nometals],
                                         "Zero_Valent":[nometals],
                                         "noint_flag":[nometals],
                                         "low_prob_1":[nometals],
                                         "low_prob_2":[nometals],
                                         "low_prob_3":[nometals],
                                         "low_prob_multi":[nometals],
                                         "high_count":[nometals],
                                         "low_count":[nometals],})
          continue  
        #print ('identifying binding sites')
        ligand_sites = get_ligand_sites(metalsites, uniquesites)
        binding_sites = get_binding_sites(metalsites, uniquesites)
        #Now get the localized oxidation state contribution of each atom
        #need delocalized bond contributions
        dVBO = delocalisedLBO(mol)
        #then need aromatic bond contributions
        rVBO = ringVBOs(mol)
        #finally combine delocal/aromatic bond conrtibutions with localized bonding
        AON = iVBS_Oxidation_Contrib(uniquesites, rVBO, dVBO)
        #Previous only assigns an oxidation contribution to unique images of atoms,
        #also need to assign these values to redundant sites:
        rAON = redundantAON(AON, mol)
        # Split the MOF into binding domains and assign oxidation state
        # contributions to binding sites
        #print ('getting binding domains')
        binding_sphere = binding_domain(
            binding_sites, rAON, mol, usitedict)
        #print ('partitioning charge')
        bindingAON = binding_contrib(
            binding_sphere, binding_sites, rAON)
        #get metal connection network
        connected_metals = get_metal_networks(ligand_sites, binding_sphere, bindingAON)
        #get metal effective coordination number:
        mCN = getCN(ligand_sites)
        
        ONEC_inout = {}
        impossible_valence = {}
        unknown_valence = {}
        zero_valence = {}
        high_count = {}
        low_count = {}
        noint_flag = {}
        noint_balance = {}
        noint_outer = {}
        low_prob_1 = {}
        low_prob_2 = {}
        low_prob_3 = {}
        low_prob_multi = {}
        low_prob_types = set()
        pure_connect_error = False
        
        #this block assigns the oxidation state and electron count only considering binding domains
        ONEC_inner = {}
        for metal in ligand_sites:
            oxidation_state = 0
            valence = valence_e(metal)
            electron_count = valence
            impossible_valence[metal] = "GOOD"
            unknown_valence[metal] = "GOOD"
            zero_valence[metal] = "GOOD"
            high_count[metal] = "GOOD"
            low_count[metal] = "GOOD"
            noint_flag[metal] = "GOOD"
            low_prob_1[metal] = "GOOD"
            low_prob_2[metal] = "GOOD"
            low_prob_3[metal] = "GOOD"
            low_prob_multi[metal] = "GOOD"
            pure_connect_error = False
            global_error = False
            
            for ligand in ligand_sites[metal]:
                LBO = bindingAON[ligand]
                Nbridge = bridging(ligand)
                Ox = LBO/Nbridge
                oxidation_state += Ox
                if Ox >= 2:
                    mCN[metal] += 1
                if Ox >= 3:
                    mCN[metal] += 1
        
            electron_count = valence + (2 * mCN[metal]) - oxidation_state
            ONEC_inner[metal] = [oxidation_state, electron_count]
            
            
        #redistribute ONEC within metal networks based on ionization energy
        noint_balance = distribute_ONEC(ONEC_inner, connected_metals, KnownIE, ONProb, HighestKnownON, mCN)
        
        #determine and distribute outer sphere charges
        OSD = outer_sphere_domain (uniquesites, binding_sphere)
        OSC = outer_sphere_contrib(OSD, rAON)
        noint_outer = distribute_OuterSphere(noint_balance, OSC, KnownIE, ONProb, HighestKnownON, mCN)
        ONEC_inout = distribute_OuterSphere(ONEC_inner, OSC, KnownIE, ONProb, HighestKnownON, mCN)
        globaldis = global_charge_distribution(ONEC_inout, KnownIE, ONProb, HighestKnownON, mCN)
            
        #flagging routine starts here.
        #start with easy flags first (impossible, unknown, zero, single improbable)
        #easy flags are based on a single criteria from a single metal
        for metal in noint_outer:
            # if the redistributed oxidation state is greater than 
            # metal valence, flag as impossible valence
            valence = valence_e(metal)
            if noint_outer[metal][0] > valence:
                impossible_valence[metal] = "BAD"   
            # if the distributed oxidation state is zero
            # flag as zero valent
            if noint_outer[metal][0] == 0:
                zero_valence[metal] = "BAD"
            # if the oxidation state after distribution is not within 0.5
            # of a known oxidation state, flag as unknown
            if not any([math.isclose(
                                    float(noint_outer[metal][0]),
                                        float(i),
                                        abs_tol = 0.5
                                        ) for i in KnownON[(metal.atomic_symbol)]
                                        ]):
                unknown_valence[metal] = "BAD"
            # if the distributed oxidation state is unlikely
            # AND a whole number, flag as low-probability
            if any([
                (math.isclose((noint_outer[metal][0] % 1), 0, abs_tol=0.0001)),
                (math.isclose((noint_outer[metal][0] % 1), 1, abs_tol=0.0001))
                ]):
                try:
                    prob = float(ONProb[metal.atomic_symbol][round(noint_outer[metal][0])])
                except IndexError:
                    prob = 0
                #metal in a less than 1% occurence oxidation state?
                #LIKELY_BAD, don't flag as hard
                if prob < 0.01:
                    low_prob_1[metal] = "LIKELY_BAD"
                    low_prob_types.add(metal.atomic_symbol)
                #less than 0.1% occurence? Flag hard.
                if prob < 0.001:
                    low_prob_2[metal] = "BAD"
                #less than 0.01%? HARD FLAG.
                if prob < 0.0001:
                    low_prob_3[metal] = "BAD"
            # for non f-block, metal electron count above 20 
            # is suspicious, flag for inspection
            if (noint_outer[metal][1] > 20) and (not mendeleev.element(metal.atomic_symbol).block == 'f'):
                high_count[metal] = "INSPECT"
            # for f-block, metal electron count above 32 
            # is suspicious, flag for inspection
            if (noint_outer[metal][1] > 32) and (mendeleev.element(metal.atomic_symbol).block == 'f'):
                high_count[metal] = "INSPECT"
            # metal electron count below 14
            # is always suspicious, flag for inspection    
            if noint_outer[metal][1] < 14:
                low_count[metal] = "INSPECT"
        if len(low_prob_types) > 1:
            for metal in noint_outer:
                low_prob_multi[metal] = "BAD"
        
        #more complicated flags here
        for metal in noint_outer:
            # if the oxidation state before distribution is zero 
            # and stays below 1 after distribution,
            #flag as zero valent and non-integer
            if (ONEC_inner[metal][0] == 0) and (noint_outer[metal][0] < 1):
                zero_valence[metal] = "BAD"
                noint_flag[metal] = "BAD"
            # non-integer flags:
            if not any([
                (math.isclose((noint_outer[metal][0] % 1), 0, abs_tol=0.0001)),
                (math.isclose((noint_outer[metal][0] % 1), 1, abs_tol=0.0001))
                ]):
                # if non-integer after distribution and 
                # impossible before distribution,
                # flag non-integer as bad
                if all([(ONEC_inner[metal][0] > valence),
                        (ONEC_inout[metal][0] > valence)]):
                    noint_flag[metal] = "BAD"
                # if non-integer after distribution and
                # zero before averaging,
                # flag as bad
                if all([(ONEC_inner[metal][0] > valence),
                        (ONEC_inout[metal][0] > valence)]):
                    zero_valence[metal] = "BAD"
                # if non-integer after averaging and
                # unknown before averaging,
                # flag as bad
                if not any([math.isclose(
                                    float(ONEC_inner[metal][0]),
                                        float(i),
                                        abs_tol = 0.5
                                        ) for i in KnownON[(metal.atomic_symbol)]
                                        ]):
                    impossible_valence[metal] = "BAD"
                if any([
                    (zero_valence[metal] == "BAD"),
                    (unknown_valence[metal] == "BAD"),
                    (impossible_valence[metal] == "BAD")
                ]):
                    noint_flag[metal] = "BAD"
                else:
                    noint_flag[metal] = "LIKELY_BAD"

        for metal in ONEC_inner:
            valence = valence_e(metal)
            if ONEC_inout[metal][0] > valence:
                pure_connect_error = True
            if ONEC_inout[metal][0] == 0:
                pure_connect_error = True
            if not any([math.isclose(
                                    float(ONEC_inout[metal][0]),
                                    float(i),
                                    abs_tol = 0.5
                                    ) for i in KnownON[(metal.atomic_symbol)]
                                    ]):
                pure_connect_error = True
            if not any([
                (math.isclose((ONEC_inout[metal][0] % 1), 0, abs_tol=0.0001)),
                (math.isclose((ONEC_inout[metal][0] % 1), 1, abs_tol=0.0001))
                ]):
                pure_connect_error = True
            if any([
                (math.isclose((ONEC_inout[metal][0] % 1), 0, abs_tol=0.0001)),
                (math.isclose((ONEC_inout[metal][0] % 1), 1, abs_tol=0.0001))
                ]):
                try:
                    prob = float(ONProb[metal.atomic_symbol][round(ONEC_inout[metal][0])])
                except IndexError:
                    prob = 0
                if prob < 0.01:
                    pure_connect_error = True
                #less than 0.1% occurence? Flag hard.
                if prob < 0.001:
                    pure_connect_error = True
                #less than 0.01%? HARD FLAG.
                if prob < 0.0001:
                    pure_connect_error = True
        if pure_connect_error == False:
            for metal in ONEC_inner:
                impossible_valence[metal] = "GOOD"
                unknown_valence[metal] = "GOOD"
                zero_valence[metal] = "GOOD"
                noint_flag[metal] = "GOOD"
                low_prob_1[metal] = "GOOD"
                low_prob_2[metal] = "GOOD"
                low_prob_3[metal] = "GOOD"
                low_prob_multi[metal] = "GOOD"
                
        #check against global charge distribution
        for metal in ONEC_inner:
            valence = valence_e(metal)
            if globaldis[metal][0] > valence:
                global_error = True
            if globaldis[metal][0] == 0:
                global_error = True
            if not any([math.isclose(
                                    float(globaldis[metal][0]),
                                    float(i),
                                    abs_tol = 0.5
                                    ) for i in KnownON[(metal.atomic_symbol)]
                                    ]):
                global_error = True
            if not any([
                (math.isclose((globaldis[metal][0] % 1), 0, abs_tol=0.0001)),
                (math.isclose((globaldis[metal][0] % 1), 1, abs_tol=0.0001))
                ]):
                global_error = True
            if any([
                (math.isclose((globaldis[metal][0] % 1), 0, abs_tol=0.0001)),
                (math.isclose((globaldis[metal][0] % 1), 1, abs_tol=0.0001))
                ]):
                try:
                    prob = float(ONProb[metal.atomic_symbol][round(globaldis[metal][0])])
                except IndexError:
                    prob = 0
                if prob < 0.01:
                    global_error = True
                #less than 0.1% occurence? Flag hard.
                if prob < 0.001:
                    global_error = True
                #less than 0.01%? HARD FLAG.
                if prob < 0.0001:
                    global_error = True
        if global_error == False:
            for metal in ONEC_inner:
                impossible_valence[metal] = "GOOD"
                unknown_valence[metal] = "GOOD"
                zero_valence[metal] = "GOOD"
                noint_flag[metal] = "GOOD"
                low_prob_1[metal] = "GOOD"
                low_prob_2[metal] = "GOOD"
                low_prob_3[metal] = "GOOD"
                low_prob_multi[metal] = "GOOD"
            
        for metal in noint_outer:
            entrydata = pd.DataFrame({"CIF":[file],
                                         "Metal":[metal.label], 
                                         "ON_coordination_ONLY":[ONEC_inner[metal][0]],
                                         "EC_coordination_ONLY":[ONEC_inner[metal][1]],
                                         "ON_coordination+Outer_Sphere":[ONEC_inout[metal][0]],
                                         "EC_coordination+Outer_Sphere":[ONEC_inout[metal][1]],
                                         "ON_network_redistributed":[noint_balance[metal][0]],
                                         "EC_network_redistributed":[noint_balance[metal][1]],
                                         "ON_global":[globaldis[metal][0]],
                                         "EC_global":[globaldis[metal][1]],
                                         "ON_network+Outer_Sphere":[noint_outer[metal][0]],
                                         "EC_network+Outer_Sphere":[noint_outer[metal][1]],
                                         "Impossible":[impossible_valence[metal]],
                                         "Unknown":[unknown_valence[metal]],
                                         "Zero_Valent":[zero_valence[metal]],
                                         "noint_flag":[noint_flag[metal]],
                                         "low_prob_1":[low_prob_1[metal]],
                                         "low_prob_2":[low_prob_2[metal]],
                                         "low_prob_3":[low_prob_3[metal]],
                                         "low_prob_multi":[low_prob_multi[metal]],
                                         "high_count":[high_count[metal]],
                                         "low_count":[low_count[metal]],})
            alldata = pd.concat([alldata,entrydata])
        alldata.to_csv(os.path.join(cwd,'OxStatesOutput.csv'))
        if not len(ZeroVectorlist) == 0:
          failmofs = open('ZeroVectorError.txt', 'w')
          failmofs.writelines(ZeroVectorlist)

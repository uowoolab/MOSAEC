# Automatic oxidation state assigner v 1.4
# developed with Python 3.7.0
# dependencies: pandas, the CCDC python API and Mendeleev libraries
# run in same folder as .cif files of interest
# folder must also contain a KnownON.csv

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
        cif.assign_bonds()
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
        for i, line in enumerate(lines):
            if '_atom_site_fract_z' in line:
                start = i
        counting = {}
        cutoff = {}
        to_remove = []
        for i in range(start, len(lines)):
            if '.' in lines[i]:
                split = lines[i].split()
                if not split[1] in counting:
                    counting[split[1]] = 1
                elif split[1] in counting:
                    counting[split[1]] +=1
                newlabel = f'{split[1]}{counting[split[1]]}'
                if newlabel in elamnt:
                    cutoff[split[1]] = counting[split[1]]
                lines[i] = lines[i].replace(split[0], newlabel)
                if split[1] in cutoff:
                    if counting[split[1]] > cutoff[split[1]]:
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

def get_unique_sites(mole, asymmole):
    #blank list for unique sites
    uniquesites = []
    asymmcoords = []
    molecoords = []
    duplicates = []
    for atom in asymmole.atoms:
        asymmcoords.append(atom.coordinates)
    for atom in mole.atoms:
        if atom.coordinates in asymmcoords:
            if not atom.coordinates in molecoords:
                uniquesites.append(atom)
                molecoords.append(atom.coordinates)
            else:
                duplicates.append(atom)
    if len(duplicates) >= 1:
        for datom in duplicates:
            for atom in uniquesites:
                if datom.coordinates == atom.coordinates:
                    if len(datom.neighbours) > len(atom.neighbours):
                        uniquesites.remove(atom)
                        uniquesites.append(datom)
    return uniquesites

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

def assign_VBS(atom, aVBO, dVBO):
    """This function will assign a Valence-Bond-Sum (VBS) to an atom. 
    Takes one CCDC atom object, list of binding sites, and the metal-free
    molecule object as inputs"""
    VBO = 0
    if atom.is_metal:
        return(0)
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

def arom_crawl (atomlist):
    """returns a list of all atoms in an aromatic system
    takes at least one aromatic bond as input"""
    for aromatom in atomlist:
        for bond in aromatom.bonds:
            if bond.bond_type == 'Aromatic':
                for member in bond.atoms:
                    if all([
                        (not member in atomlist),
                        (not member.is_metal)
                    ]):
                        atomlist.append(member)
    return(atomlist)

def noarom_VBS(atom):
    """Valence Bond Sum of an atom excluding aromatic bonds
    necessary to determine aromatic VBS"""
    VBO = 0
    if atom.is_metal:
        return(VBO)
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

def aromaticVBO(molecule):
    """returns a dict of all aromatic atoms
    and their (aromatic-only) VBS"""
    aromdict = {}
    for atom in molecule.atoms:
        if all([
            (any(
                bond.bond_type == 'Aromatic'
                for bond in atom.bonds
                )),
                (not atom in aromdict)
            ]):
            arom_atoms = arom_crawl(atomlist=[atom])
            aromcoords = set()
            for member in arom_atoms:
                aromcoords.add(member.coordinates)
            
            #kekulize
            ringmol = molecule.copy()
            for atom in ringmol.atoms:
                if atom.is_metal:
                    ringmol.remove_atom(atom)
            ringmol.kekulize()
            
            #get charge
            ringcharge = 0
            for member in ringmol.atoms:
                if member.coordinates in (aromcoords):
                    charge = iVBS_FormalCharge(member)
                    ringcharge += charge
            
            for member in arom_atoms:
                Narombonds = 0
                for abond in member.bonds:
                    if abond.bond_type == 'Aromatic':
                        Narombonds += 1
                unpaired_e = (4 - abs(4 - valence_e(member)))
                #aromatic bond order is then:
                aVBS = (
                    ((unpaired_e-noarom_VBS(member))/Narombonds)
                    + (ringcharge/(Narombonds*len(arom_atoms))))
                aromdict[member] = aVBS
    return (aromdict)

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
    unpaired_electrons = (4 - abs(4 - valence))
    #also require "available bonding electrons"
    ABE = ((VBO-valence)*2)+valence
    if VBO > valence:
        #expanded valence centres are handled by
        charge = VBO - ABE
    if VBO <= valence:
        #normal valence centres are handled by
        charge = VBO - unpaired_electrons
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
        if atom.atomic_number in range(57,71):
            valence = atom.atomic_number - 57 + 3
        elif atom.atomic_number in range(89,103):
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
    return bridge

def iVBS_Oxidation_Contrib(unique_atoms, aVBO, dVBO):
    """determines the oxidation state contribution for all
    unique atoms in a MOF. Returns a dictionary of Atom:Oxidaiton_Contribution
    pairs. Takes the unique sites in the MOF (without metal), the MOF molecule
    object (without metal) and a list of metal-bindin sites as inputs"""
    VBS = 0
    CN = 0
    valence = 0
    oxi_contrib = {}
    # for each unique atom
    for atom in unique_atoms:
        # assign valence-bond-sum
        VBS = assign_VBS(atom, aVBO, dVBO)
        #determine coordination number
        CN = get_CN(atom)
        #  determine number of valence electrons
        valence = valence_e(atom)
        # get number of unpaired electrons in the free element
        unpaired_e = (4 - abs(4 - valence))
        #Available Bonding Electrons
        ABE = ((VBS-valence)*2)+valence
        #  Expanded valences are handled by:
        if  atom.is_metal:
            oxi_contrib[atom] = 0
        elif VBS > valence:
            oxi_contrib[atom] = ABE - VBS
        # Otherwise, normal valences:
        elif VBS <= valence:
            oxi_contrib[atom] = unpaired_e - VBS
        # need to check for 3-coordinate carbocations,
        # 3-coordinate carbanions, carbenes, and heavier
        # homologues (these are not immediately detectable)
        if any([
            (atom.atomic_symbol == 'C'),
            (atom.atomic_symbol == 'Si'),
            (atom.atomic_symbol == 'Ge'),
            (atom.atomic_symbol is 'Pb')]):
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
    return (oxi_contrib)

def redundantAON(AON, molecule):
    redAON = {}
    for rsite1 in molecule.atoms:
        for usite1 in AON:
            redAON[usite1] = AON[usite1]
            if rsite1.label == usite1.label:
                redAON[rsite1] = AON[usite1]
    return(redAON)

def binding_domain(binding_sites, AON, molecule):
    """Due to geometric distortions (as expected in a MOF),
    delocalized bonds are not always correctly assigned
    This function will build 'bonding domains' to fix this 
    (see methodology section of paper)"""
    def arom_domains(site, aromlist, bondset):
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
                    return(arom_domains(site, aromlist, bondset))
        return(aromlist)
    def deloc_domains(site, AON, molecule, deloclist, bondset):
        """recursive crawler makes domains"""
        for bond in site.bonds:
            bondset.add(bond)
        for bond in bondset:
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
                        (not mbond.bond_type == 'Single'
                        for mbond in member.bonds),
                        (not mbond.bond_type == 'Aromatic'
                        for mbond in member.bonds),
                        (not mbond.bond_type == 'Delocalised'
                        for mbond in member.bonds)
                        ]))
                    ]))
                    ])
                    ):
                    deloclist.append(member)
                    for mbond in member.bonds:
                        bondset.add(mbond)
                    return(deloc_domains(site, AON, molecule, deloclist, bondset))
        return(deloclist)
    sitedomain={}
    for site in binding_sites:
        if any(sbond.bond_type == 'Aromatic' for sbond in site.bonds):
            sitedomain[site] = arom_domains(site, aromlist=[site], bondset=set())
        if not any(sbond.bond_type == 'Aromatic' for sbond in site.bonds):
            sitedomain[site] = deloc_domains(site, AON, molecule, deloclist=[site], bondset=set())
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
                                            return network_crawl(ligand_sites, binding_sphere, bindingAON, metal_networks, checked_sites, group)
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
            metal_networks.append(network_crawl(ligand_sites, binding_sphere, bindingAON, metal_networks, checked_sites=[], group=[metal]))

    network_dict = {}
    for network in metal_networks:
        for metal in network:
            network_dict[metal] = network
    return (network_dict)

def average_noints(ONEC, metal_connection):
    fONEC = {}
    for metal in ONEC:
        sumON = []
        sumEC = []
        ONsum = 0
        ECsum = 0
        for metals in metal_connection[metal]:
            sumON.append(ONEC[metals][0])
            ONsum = math.fsum(sumON)
            sumEC.append(ONEC[metals][1])
            ECsum = math.fsum(sumEC)
        avgON = ONsum/len(metal_connection[metal])
        avgEC = ECsum/len(metal_connection[metal])
        fONEC[metal] = (avgON, avgEC)
    return (fONEC)
    
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

cwd = os.getcwd()
KnownON = KnownONs()
alldata = pd.DataFrame()
for file in os.listdir(cwd):
    if (str(file).endswith('.cif') == True):
        print (f'Analyzing {file}...')
        input_cif = os.path.join(cwd,file)
        # read in the .cif, extract the underlying molecule,
        # identify the unique sites, metal sites, binding sites,
        # ect.
        #print ('reading .cif data')
        cif = readentry(input_cif)
        mol = cif.molecule
        asymmol = cif.asymmetric_unit_molecule
        #print ('identifying unique sites')
        uniquesites = get_unique_sites(mol, asymmol)
        #print ('getting metal sites')
        metalsites = get_metal_sites(uniquesites)
        #print ('identifying binding sites')
        ligand_sites = get_ligand_sites(metalsites, uniquesites)
        binding_sites = get_binding_sites(metalsites, uniquesites)
        #Now get the localized oxidation state contribution of each atom
        #need delocalized bond contributions
        dVBO = delocalisedLBO(mol)
        #then need aromatic bond contributions
        aVBO = aromaticVBO(mol)
        #finally combine delocal/aromatic bond conrtibutions with localized bonding
        AON = iVBS_Oxidation_Contrib(
            uniquesites, aVBO, dVBO)
        #Previous only assigns an oxidation contribution to unique images of atoms,
        #also need to assign these values to redundant sites:
        rAON = redundantAON(AON, mol)
        # Split the MOF into binding domains and assign oxidation state
        # contributions to binding sites
        #print ('getting binding domains')
        binding_sphere = binding_domain(
            binding_sites, rAON, mol)
        #print ('partitioning charge')
        bindingAON = binding_contrib(
            binding_sphere, binding_sites, rAON)
        #get metal connection network
        connected_metals = get_metal_networks(ligand_sites, binding_sphere, bindingAON)
        #print ('assigning oxidation state')
        
        ONEC_with_outer = {}
        impossible_valence = {}
        unknown_valence = {}
        zero_valence = {}
        high_count = {}
        low_count = {}
        noint_flag = {}
        noint_balance = {}
        noint_outer = {}
        
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
            for ligand in ligand_sites[metal]:
                LBO = bindingAON[ligand]
                hapto = hapticity(ligand, metal)
                Nbridge = bridging(ligand)
                Ox = LBO/Nbridge
                Ec = 2/Nbridge
                if LBO < 2:
                    Ec = (2 - Ox)
                if hapto == True:
                    Ec = 1/Nbridge
                oxidation_state += Ox
                electron_count += Ec
            ONEC_inner[metal] = [oxidation_state, electron_count]
 
        noint_balance = average_noints(ONEC_inner, connected_metals)
        OSD = outer_sphere_domain (uniquesites, binding_sphere)
        OSC = outer_sphere_contrib(OSD, rAON)
        OSCa = OSC

        #this block assigns oxidation state and electron count including outer sphere domain
        for metal in ONEC_inner:
            ONEC_with_outer[metal] = list(ONEC_inner[metal])
        if not OSCa == 0:
            metalgroups = equiv_metals(metalsites, rAON,  ONEC_with_outer)
            if len(metalgroups) == 1:
            # if there is only one type of metal, distribute
            # outer sphere charge equally to all metals
                for metal in ONEC_with_outer:
                    ONEC_with_outer[metal][0] += OSCa/len(ONEC_with_outer)
                    ONEC_with_outer[metal][1] -= OSCa/len(ONEC_with_outer)
                OSCa = 0
            if len(metalgroups) > 1:
                # if there is more than one type of metal, we have a few things to try
                # first, if there are non-averaging non-integer oxidation states, 
                # add charge to bring these to integer value 
                for metal in ONEC_with_outer:
                    int_gap = (1 - (ONEC_with_outer[metal][0] % 1))
                    mod_outer = (ONEC_with_outer[metal][0] % 1)
                    mod_balance = (noint_balance[metal][0] % 1)
                    if all([(abs(OSCa) >= int_gap),
                            (not math.isclose(mod_outer , 0)),
                            (not math.isclose(mod_balance, 0)),
                            (not math.isclose(mod_outer , 1)),
                            (not math.isclose(mod_balance, 1)),
                            ]):
                        if OSCa > 0:
                            OSCa -= int_gap
                            ONEC_with_outer[metal][1] -= int_gap
                            ONEC_with_outer[metal][0] += int_gap
                        if OSCa < 0:
                            OSCa += mod_outer
                            ONEC_with_outer[metal][1] += mod_outer
                            ONEC_with_outer[metal][0] -= mod_outer
                # now that all metals (ideally) have integer oxidation states,
                # (or non-integer oxidation states that balance)
                # if there are impossible valences (and there is enough outer-sphere
                # to do so), fix those first
                for metal in ONEC_with_outer:
                    valence = valence_e(metal)
                    diff = (ONEC_with_outer[metal][0] - valence)
                    if all([
                        (abs(OSCa) >= abs(diff)),
                        (ONEC_with_outer[metal][0] > valence),
                        (OSCa < 0),
                        ]):
                        OSCa += diff
                        ONEC_with_outer[metal][0] -= diff
                        ONEC_with_outer[metal][1] += diff      
                #if there are still charges to assign, split evenly
                if not math.isclose(OSCa, 0) :
                    for metal in ONEC_with_outer:
                        ONEC_with_outer[metal][0] += OSCa/len(ONEC_with_outer)
                        ONEC_with_outer[metal][1] -= OSCa/len(ONEC_with_outer)
                    OSCa = 0
            
        noint_outer = average_noints(ONEC_with_outer, connected_metals)
        for metal in ONEC_with_outer:
            #quick check for impossible or unknown valences
            if noint_outer[metal][0] >= (valence + 0.5):
                impossible_valence[metal] = "BAD"
            if noint_outer[metal][0] == 0:
                zero_valence[metal] = "BAD"
            if not any([math.isclose(
                                    float(noint_outer[metal][0]),
                                        float(i),
                                        abs_tol = 0.5
                                        ) for i in KnownON[metal.atomic_symbol]
                                        ]):
                unknown_valence[metal] = "BAD"
            if noint_outer[metal][0] < 0.9:
                unknown_valence[metal] = "BAD"
            if noint_outer[metal][1] >= 20:
                high_count[metal] = "INSPECT"
            if noint_outer[metal][1] <= 14:
                low_count[metal] = "INSPECT"
            if not any([
                        (math.isclose((noint_outer[metal][0] % 1), 0, abs_tol=0.0001)),
                        (math.isclose((noint_outer[metal][0] % 1), 1, abs_tol=0.0001))]):
                if any([
                    (zero_valence[metal] == "BAD"),
                    (unknown_valence[metal] == "BAD"),
                    (impossible_valence[metal] == "BAD")
                ]):
                    noint_flag[metal] = "BAD"
                else:
                    noint_flag[metal] = "LIKELY BAD"

            entrydata = pd.DataFrame({"CIF":[file],
                                         "Metal":[metal.label], 
                                         "ON":[ONEC_with_outer[metal][0]],
                                         "EC":[ONEC_with_outer[metal][1]],
                                         "ON_coordination_ONLY":[ONEC_inner[metal][0]],
                                         "EC_coordination_ONLY":[ONEC_inner[metal][1]],
                                         "ON_noint_average":[noint_outer[metal][0]],
                                         "EC_noint_average":[noint_outer[metal][1]],
                                         "Impossible":[impossible_valence[metal]],
                                         "Unknown":[unknown_valence[metal]],
                                         "Zero_Valent":[zero_valence[metal]],
                                         "noint_flag":[noint_flag[metal]],
                                         "high_count":[high_count[metal]],
                                         "low_count":[low_count[metal]],})
            alldata = pd.concat([alldata,entrydata])
        alldata.to_csv(os.path.join(cwd,'OxStatesOutput.csv'))
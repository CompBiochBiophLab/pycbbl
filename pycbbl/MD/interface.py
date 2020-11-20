import mdtraj as md
import numpy as np

def getInterfaceResidues(trajectory, frame, threshold_sasa=0.2, exclude=[],
                         return_hsasa=False, return_phsasa=False):
    """
    Calculate the interface residues between the different chains in the model.
    The function returns a dictionary having the chain indexes of a particular interface
    as keys and the residues belonging to that interface as its values.
    The function only works with one frame for easy reporting, but it could be modified
    if needed to include all frames in the analysis.

    Parameters
    ----------
    trajectory : mdtraj.Trajectory
        Multichain trajetory to analyse
    frame : int
        Frame to calculate residues in interface
    threshold_sasa : float
        Percentage of hidden SASA to select residues
    exclude : list (int)
        List of chain indexes to exclude from the analysis
    return_hsasa : bool
        Return hidden sasa per-residue instead of interface residues
    Returns
    -------
    interface_residues : dict
        Dictionary containing the residues in all the calculated interfaces.
    """

    # Save original trajetory for referencing residues
    topology = trajectory.topology
    # Slice trajectory to extract specific frame and create a new copy of the object
    trajectory = trajectory.slice(frame)
    # Get trajectory chains as dictionary accesible by chain index
    chains = {c.index:c for c in trajectory.topology.chains}
    if len(chains) <= 1:
        raise ValueError('Trajectory topology does not contain more than one chain.\
The input trajectory should be multi chain trajectory object.')

    if isinstance(exclude, int):
        exclude = list(exclude)
    elif not isinstance(exclude, list):
        raise ValueError('Excluded chains should and integer or a list of integers')

    # Store atom indexes by chain
    chain_atoms = {}
    for c in chains:
        chain_atoms[c] = [a.index for a in chains[c].atoms]

    # Create dictionary to store interface residues
    interface_residues = {}
    # Create dictionary to store residues into dictionary
    residues = {}

    # Check which chains will be compared
    for c in chains:
        for d in chains:
            if c != d and c < d:
                if c not in exclude and d not in exclude:
                    # Save indexes into interface residues dictionary
                    interface_residues[(c, d)] = []
                    # Save original topology residues into dictionary
                    residues[(c, d)] = []
                    for cx in topology.chains:
                        if cx.index in (c, d):
                            for r in cx.residues:
                                residues[(c, d)].append(r)

    hidden_sasa = {}
    phidden_sasa = {}

    # Calculate solvent accesible areas for complex and single chains
    for c, d in interface_residues:
        hidden_sasa[(c,d)] = {}
        phidden_sasa[(c,d)] = {}
        # Create trajectories containing all chain atoms and single chain atoms
        complex_atoms = chain_atoms[c] + chain_atoms[d]
        t = trajectory.atom_slice(complex_atoms)
        tc = trajectory.atom_slice(chain_atoms[c])
        td = trajectory.atom_slice(chain_atoms[d])

        # Renumber residues to lift restriction of the shrake_rupley method
        renumberContiguousResidues(t)
        renumberContiguousResidues(tc)
        renumberContiguousResidues(td)

        # Calculate complex and single chain SASAs
        complex_sasa = md.shrake_rupley(t, mode='residue')
        c_sasa = md.shrake_rupley(tc, mode='residue')
        d_sasa = md.shrake_rupley(td, mode='residue')

        # Put together single chain SASA arrays to compare with complex SASA
        unbound_sasa = np.concatenate((c_sasa, d_sasa), axis=1)

        # Compare complex and unbound SASA and save residues with hidden sasa larger then threshold
        hsasa = unbound_sasa - complex_sasa
        hidden_percentage = hsasa / unbound_sasa

        # Get residues with hidden SASA larger than threshold_sasa (default 20%)
        for i,v in enumerate(hidden_percentage[0]):
            # Sum hidden sasa for each residue at each interface
            hidden_sasa[(c,d)][residues[c,d][i]] = hsasa[0][i]
            phidden_sasa[(c,d)][residues[c,d][i]] = hidden_percentage[0][i]
            if v > threshold_sasa:
                interface_residues[(c, d)].append(residues[c,d][i])

    if return_hsasa:
        return hidden_sasa
    elif return_phsasa:
        return phidden_sasa
    return interface_residues

def calculateSasa(trajectory, chain_indexes=[], mode='residue'):
    """
    Caluculate the SASA for the full or a subset of chains in a trajectory.

    Parameters
    ----------
    chain_indexes : list
        Optional list of chain indexes to consider in the SASA calculation.

    Returns
    -------
    sasa : numpy.ndarray
        SASA values for the residues in the trajectory chains.
    """

    # Select all chain indexes if none were given.
    if chain_indexes == []:
        for c in trajectory.topology.chains:
            chain_indexes.append(c.index)

    # Select atoms in selected chains
    chains_atoms = []
    for c in trajectory.topology.chains:
        if c.index in chain_indexes:
            for a in c.atoms:
                chains_atoms.append(a.index)

    # Create trajectory with subset of selected atoms
    t = trajectory.atom_slice(chains_atoms)
    renumberContiguousResidues(t)
    sasa = md.shrake_rupley(t, mode='residue')

    return sasa

def renumberContiguousResidues(trajectory):

    for i,r in enumerate(trajectory.topology.residues):
        r.index = i

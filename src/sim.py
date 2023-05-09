

##
##@Author: Greg Pintilie
##@Date:   14-Apr-2023
##@Email:  gregdp@gmail.com
##@Last modified by:   gregdp
##@Last modified time: 14-Apr-2023
##@License: Free for non-commercial use (see license.pdf)
##@Copyright: Copyright 2023 Grigore Pintilie
##



openmm_forcefield_parameters = ['amber14-all.xml', 'amber14/tip3p.xml']
write_logs = True

from .sim_mm_pull import PullAtomMode, Logger

class MolSim :

    '''
    Molecular dynamics simulation with OpenMM
    '''

    def __init__(self, atoms, dmap=None, force_constant = 10000.0,
                 cutoff = 10.0, temperature = 100.0,
                 steps = 50, tolerance = 0.001):

        self._log = Logger('/Users/greg/Desktop/sim.log' if write_logs else None)
        self._log ( "Initialized MolSim" )
        self._mol = atoms[0].structure
        self._atoms = atoms
        self._minimized = False
        self._dmap = dmap

        initialize_openmm()

        # OpenMM objects
        self._topology = None
        self._system = None
        self._force = None	# CustomExternalForce pulling force
        self._platform = None
        self._simulation = None
        self._sim_forces = None # Current forces on atoms
        self._integrator = None

        from time import time
        self._sim_t = time()


        # OpenMM simulation parameters
        global openmm_forcefield_parameters
        self._forcefields = openmm_forcefield_parameters
        self._sim_steps = steps		# Simulation steps between mouse position updates
        self.force_constant = force_constant
        self.dforce_constant = force_constant
        self._nonbonded_cutoff = cutoff
        from openmm import unit
        self._temperature = temperature * unit.kelvin
        self._integrator_tolerance = tolerance
        self._constraint_tolerance = tolerance
        self._friction = 1.0/unit.picoseconds	# Coupling to heat bath
        self._platform_name = 'CPU'
        #self._platform_name = 'OpenCL' # Works on Mac
        #self._platform_name = 'CUDA'	# This is 3x faster but requires env DYLD_LIBRARY_PATH=/usr/local/cuda/lib Chimera.app/Contents/MacOS/ChimeraX so paths to cuda libraries are found.
        self._max_allowable_force = 5e6 # kJ/mol/nm


        # OpenMM particle data
        self._particle_positions = None		# Positions for all atoms, numpy array, Angstroms
        self._particle_force_index = {}		# Maps particle number to force index for tugged atoms
        self._particle_masses = None		# Original particle masses

        self.pull_atom_i = None
        self.pull_atom_pt = None
        self.pull_atom = False

        self._create_openmm_system()



    def start ( self ) :
        print ( " - make simulation" )
        if self._simulation != None :
            print ( " - ?" )
        else :
            self._make_simulation()
            print ( " - made simulation" )


    def step ( self, num_steps = None ) :
        self._simulate ( steps=num_steps )


    def SetTemp ( self, temp ) :
        if self._integrator :
            print ( " - setting temp %.3f" % temp )
            from openmm import unit
            self._temperature = temp * unit.kelvin
            self._integrator.setTemperature ( self._temperature )


    def _make_simulation ( self ) :
        self._log ( " _ making simulation" )
        import openmm as mm
        self._integrator = mm.VariableLangevinIntegrator(self._temperature, self._friction, self._integrator_tolerance)
        self._integrator.setConstraintTolerance(self._constraint_tolerance)

        # Make a new simulation.
        from openmm import app
        s = app.Simulation(self._topology, self._system, self._integrator, self._platform)
        self._simulation = s


    def _max_force ( self ):
        c = self._simulation.context
        from openmm.unit import kilojoule_per_mole, nanometer
        self._sim_forces = c.getState(getForces = True).getForces(asNumpy = True)/(kilojoule_per_mole/nanometer)
        forcesx = self._sim_forces[:,0]
        forcesy = self._sim_forces[:,1]
        forcesz = self._sim_forces[:,2]
        import numpy
        magnitudes = numpy.sqrt(forcesx*forcesx + forcesy*forcesy + forcesz*forcesz)
        return max(magnitudes)

    def _get_energy ( self ) :
        c = self._simulation.context
        self._k_energy = c.getState(getEnergy=True).getKineticEnergy()
        self._p_energy = c.getState(getEnergy=True).getPotentialEnergy()

        from openmm.unit import kilojoule_per_mole
        self._k_energy = self._k_energy.value_in_unit(kilojoule_per_mole)
        self._p_energy = self._p_energy.value_in_unit(kilojoule_per_mole)

        return self._p_energy, self._k_energy



    def pull_atom_to_pt (self, atom, atom_i, xyz):

        '''xyz is in scene coordinates, one position for each tugged particle.'''
        if 0 :
            self._log( " tug_to_positions - %d " % len(xyz)  )
            for i in range ( len(xyz) ) :
                self._log ( " [%d] : %.3f, %.3f %.3f" % (i, xyz[i][0], xyz[i][1], xyz[i][2]) )

        txyz = self._mol.scene_position.inverse() * xyz
        #txyz *= 0.1	# convert position in Angstroms to nanometers
        v = txyz - atom.coord

        self._log ( " - old - force for atom %s %d -> %.3f,%.3f,%.3f" % (atom.name, atom_i, v[0], v[1], v[2]) )
        f = self._force
        k = self.force_constant
        f.setParticleParameters(atom_i, atom_i, (k,v[0],v[1],v[2]))
        f.updateParametersInContext(self._simulation.context)


        #if not self._minimized:
        #    self._minimize()
        #self._simulate()

    def reset_forces ( self ) :

        f = self._force
        for ai, at in enumerate ( self.atoms ) :
            f.setParticleParameters(ai, ai, (0,0,0,0))
        f.updateParametersInContext(self._simulation.context)



    def _set_map_forces ( self ) :

        from time import time, sleep

        max_map_g = 0.0
        kp = self.force_constant
        kd = self.dforce_constant
        f = self._force
        max_k = self._max_allowable_force

        if self.pull_atom_i != None and self.pull_atom != False :
            atom, atom_i = self.atoms[self.pull_atom_i], self.pull_atom_i
            txyz = self._mol.scene_position.inverse() * self.pull_atom_pt
            v = txyz - atom.coord
            #self._log ( " - pull force for atom %s %d -> %.3f,%.3f,%.3f - kp %f" % (atom.name, atom_i, v[0], v[1], v[2], kp ) )
            f.setParticleParameters(atom_i, atom_i, (kp,v[0],v[1],v[2]))

        if self._dmap != None :
            #vals = self._dmap.interpolated_values ( self._hatoms.scene_coords )
            dmapM, molM = self._dmap.scene_position, self._mol.scene_position
            atPts = self._hatoms.scene_coords
            atPts = dmapM.inverse() * atPts
            grads = self._dmap.interpolated_gradients ( atPts ) # in mol coords with xform faster?
            grads = (molM.inverse() * dmapM).transform_vectors (grads)
            #t2 = time()

            gx, gy, gz = grads[:,0], grads[:,1], grads[:,2]
            from numpy import sqrt
            magnitudes = sqrt(gx*gx + gy*gy + gz*gz)
            max_map_g = max(magnitudes)

            if kd * max_map_g > self._max_allowable_force * 0.5 :
                kd = self._max_allowable_force / max_map_g * 0.5
                self._log ( " - reduced kd %f -> %f" % (self.dforce_constant, kd) )

            if self.pull_atom_i == None or self.pull_atom == False :
                for i, g in enumerate (grads) :
                    atom_i = self._hiatoms[i]
                    f.setParticleParameters (atom_i, atom_i, (kd,g[0],g[1],g[2]))

                #self._log ( " - set %d forces (max %f / kd %f (%f), max_k %f) atoms" % (len(grads), maxg, kd, maxg*kd, max_k) )
            elif 0 :
                for i, g in enumerate ( grads ) :
                    atom_i = self._hiatoms[i]
                    atom = self._hatoms[i]
                    if atom_i == self.pull_atom_i :
                        txyz = self._mol.scene_position.inverse() * self.pull_atom_pt
                        v = txyz - atom.coord
                        #self._log ( " - pull force for atom %s %d -> %.3f,%.3f,%.3f - kp %f" % (atom.name, atom_i, v[0], v[1], v[2], kp ) )
                        f.setParticleParameters(atom_i, atom_i, (kp,v[0],v[1],v[2]))
                    else :
                        f.setParticleParameters (atom_i, atom_i, (kd,g[0],g[1],g[2]))
            elif 1 :
                gi = 0
                for atom_i, atom in enumerate ( self.atoms ) :
                    if atom_i == self.pull_atom_i :
                        txyz = self._mol.scene_position.inverse() * self.pull_atom_pt
                        v = txyz - atom.coord
                        #self._log ( " - pull force for atom %s %d -> %.3f,%.3f,%.3f - kp %f" % (atom.name, atom_i, v[0], v[1], v[2], kp ) )
                        f.setParticleParameters(atom_i, atom_i, (kp,v[0],v[1],v[2]))
                    else :
                        if atom.element.number != 1 :
                            g = grads[gi]; gi += 1
                            f.setParticleParameters (atom_i, atom_i, (kd,g[0],g[1],g[2]))


        f.updateParametersInContext(self._simulation.context)

        return max_map_g



    def _simulate (self, steps = None) :

        from time import time

        t0 = time ()
        self._simulation.context.setPositions(0.1*self.atoms.coords)	# Nanometers
        self._simulation.context.setVelocitiesToTemperature(self._temperature)
        max_map_g = self._set_map_forces ()
        t1 = time()

        sim_steps = self._sim_steps if steps is None else steps

        t2 = time()
        sim_steps = 1
        self._simulation.step(sim_steps)
        t3 = time()

        if 1 :
            max_force = self._max_force()
            if max_force > self._max_allowable_force :
                self._log ( " - max force: %g / %g - minimizing..." % (max_force, self._max_allowable_force) )
                self.minimize()

        t4 = time()
        self._update_atom_coords()
        t5 = time()

        pE, kE = self._get_energy ()

        at_t = time()
        if at_t - self._sim_t > 2.0 :
            self._log ( ' - %d steps in %.4f/%.4f/%.4f seconds - max force %.3g/%.3g - energy p:%.2f k:%.2f - fk %g dk %g - maxg %g' % (sim_steps, t1-t0, t3-t2, t5-t4, max_force, self._max_allowable_force, pE, kE, self.force_constant, self.dforce_constant, max_map_g))
            self._sim_t = at_t



    def minimize ( self ) :

        if self._system == None :
            self._mol.session.logger.error ( "Trying to minimize but system has not been created" )
            return

        self._make_simulation()

        self.reset_forces()
        self._simulation.context.setPositions(0.1*self.atoms.coords)	# Nanometers
        #self._simulation.context.setVelocitiesToTemperature(self._temperature)
        from openmm import unit
        self._simulation.context.setVelocitiesToTemperature (0.0 * unit.kelvin)

        #max_force = self._max_force()
        #last_force = max_force

        pot0, kin0 = self._get_energy ()
        self._log ( "Minimizing..." )
        self._log ( " - potential E: %f, kinetic E: %f" % (pot0, kin0) )

        from time import time
        t0, t_at = time(), time()
        sim_steps = 10
        for ati in range ( 10000 ) :

            try :
                self._simulation.minimizeEnergy(maxIterations = sim_steps)
            except Exception as e1:
                self._log ( "Minimize exception: %s" % str(e1) )
                print ( "Minimize exception: %s" % str(e1) )

            pot, kin = self._get_energy ()
            tot = pot + kin
            diff = (pot+kin) - (pot0+kin0)

            t = time()
            if t - t_at > 2.0 :
                self._log ( " - %d - potential E: %f, kinetic E: %f, tot E: %f, change: %f" % (ati+1, pot, kin, tot, diff) )
                t_at = t
                self._update_atom_coords()

            if t - t0 > 10.0 :
                self._log ( " - %d - potential E: %f, kinetic E: %f, tot E: %f, change: %f" % (ati+1, pot, kin, tot, diff) )
                self._log ( " - stopping after 10s" )
                self._update_atom_coords()
                break

            if abs(diff) <= 1.0 :
                self._log ( " - %d - potential E: %f, kinetic E: %f, tot E: %f, change: %f" % (ati+1, pot, kin, tot, diff) )
                self._log ( " - stopping" )
                break

            pot0, kin0 = pot, kin




    def _minimize(self, steps = None):
        self._simulation.context.setPositions(0.1*self.atoms.coords)	# Nanometers
        self._simulation.context.setVelocitiesToTemperature(self._temperature)
        min_steps = self._sim_steps if steps is None else steps
        self._simulation.minimizeEnergy(maxIterations = min_steps)
        self._update_atom_coords()
        self._minimized = True

    def mobile_atoms(self, atoms):
        '''
        Fix positions of some particles.  Must be called before creating OpenMM simulation otherwise
        it has no effect.

        This works by an OpenMM convention that zero mass particles do not move.
        But the OpenMM docs says contraints to zero mass particles don't work.
        This means bond length constraints cannot be used to allow longer integration
        time steps. For reasons I do not understand, OpenMM it will work.
        '''
        np = len(self.atoms)
        m = self._particle_masses
        system = self._system
        if m is None:
            self._particle_masses = m = [system.getParticleMass(i) for i in range(np)]
        mi = set(self.atoms.indices(atoms))
        freeze_mass = 0
        for i in range(np):
            mass = m[i] if i in mi else freeze_mass
            system.setParticleMass(i, mass)


    def _update_atom_coords(self):
        c = self._simulation.context
        state = c.getState(getPositions = True)
        from openmm import unit
        pos = state.getPositions().value_in_unit(unit.angstrom)
        from numpy import array, float64
        xyz = array(pos, float64)
        self.atoms.coords = xyz


    def _create_openmm_system(self):

        #self._log('In create_openmm_system ')

        from chimerax.core.commands import run

        if 0 :
            mid = "%d" % self._mol.id[0]
            for i in range ( 1, len(self._mol.id) ) :
                mid += ".%d" % self._mol.id[i]
            cmd = "addh #%s" % (mid)
            print ( " - running addh: %s" % cmd )

            run(self._mol.session, cmd  )

        print ( " - creating system..." )

        system = None
        last_pres = None
        while 1 :

            self._log ( "Creating system..." )

            # OpenMM requires the atoms to be sorted by residue number
            satoms = list(self._atoms)
            satoms.sort(key = lambda a: (a.residue.chain_id, a.residue.number))
            from chimerax.atomic import Atoms, Bonds
            self.atoms = Atoms(satoms)

            atomsMap = {}
            for at in self.atoms :
                atomsMap[at] = 1

            bonds = [] # self._mol.bonds
            for b in self._mol.bonds :
                if b.atoms[0] in atomsMap and b.atoms[1] in atomsMap :
                    bonds.append ( b )

            self.bonds = Bonds ( bonds )

            self._topology = openmm_topology(self.atoms, self.bonds)
            from openmm import app
            forcefield = app.ForceField(*self._forcefields)
    #        self._add_hydrogens(pdb, forcefield)
            self._add_ligands_to_forcefield(self._mol.residues, forcefield)
            system, prob_res = self._create_system(forcefield)

            break

            if system == None and prob_res != None :
                if last_pres != None and prob_res.chain_id == last_pres.chain_id and prob_res.number == last_pres.number :
                    print ( " - alread tried to fix, didn't work" )
                    self._log ( " - alread tried to fix, didn't work" )
                else :
                    res = prob_res
                    cmd = None
                    if self._dmap == None :
                        cmd = "swap /%s:%d %s criteria c" % (res.chain_id, res.number, res.name)
                    else :
                        did = "%d" % self._dmap.id[0]
                        for i in range ( 1, len(self._dmap.id) ) :
                            did += ".%d" % self._dmap.id[i]
                        cmd = "swap /%s:%d %s criteria d density #%s" % (res.chain_id, res.number, res.name, did)
                    print ( " = running %s" % cmd)
                    self._log ( " = running %s" % cmd )
                    self._log ( cmd )
                    try :
                        run ( self._mol.session, cmd  )
                    except Exception as e1 :
                        print ( " - couldn't swap: %s" % str(e1) )
                        self._log ( " - couldn't swap: %s" % str(e1) )
                        break
                    last_pres = res
                    continue

            break

        if system == None :
            raise ForceFieldError('Missing atoms or parameterization needed by force field.\n' +
                                  'All heavy atoms and hydrogens with standard names are required.\n' )


        self._system = system

#        path = '/Users/goddard/ucsf/amber/sustiva_tutorial/1FKO_sus.prmtop'
#        crd_path = path[:-7] + '.rst7'
#        path = '/Users/goddard/ucsf/amber/gfp_tutorial/gfp.parm7'
#        crd_path = '/Users/goddard/ucsf/amber/gfp_tutorial/min1.rst7'
#         self._system_from_prmtop(path, crd_path)

        import openmm as mm
        platform = mm.Platform.getPlatformByName(self._platform_name)
        self._platform = platform
        #print ( " -- platform: ", platform )
        self._log ( " -- platform: %s" % platform )

        # http://docs.openmm.org/7.1.0/api-python/generated/simtk.openmm.openmm.CustomExternalForce.html
        #e = 'k*((x-x0)^2+(y-y0)^2+(z-z0)^2)'

        # https://openmm.github.io/openmm-cookbook/latest/notebooks/cookbook/Applying%20a%20Constant%20External%20Force.html
        e = 'k*(-fx*x-fy*y-fz*z)'

        self._force = force = mm.CustomExternalForce(e)

        for name in ('k', 'fx', 'fy', 'fz'):
            force.addPerParticleParameter(name)

        hatoms, self._hiatoms = [], []

        for ai, at in enumerate ( self.atoms ) :
            ri = force.addParticle(ai, (0,0,0,0))
            if ai != ri :
                self._log ( " - returned index different %d -> %d " % (ai, ri) )
                #self._mol.session.logger.error ( "unexpected result in start sim, probably will not run ok" )
                #return
            if at.element.name != "H" :
                hatoms.append ( at )
                self._hiatoms.append ( ai )

            #force.setParticleParameters(ai, ai, (0,0,0,0))

        from chimerax.atomic import Atoms
        self._hatoms = Atoms(hatoms)

        self._log ( " - added %d atoms as particles to force, %d heavy" % (len(self.atoms),len(self._hatoms)) )

        system.addForce ( self._force )
        #self._force.updateParametersInContext(self._simulation.context)

        print ( " - system created!" )
        self._log ( " - system created!" )



    def _add_ligands_to_forcefield(self, residues, forcefield):

        known_ligands = set(forcefield._templates.keys())
        known_ligands.add('HIS')	# Seems OpenMM handles this specially.
        from numpy import unique
        rnames = [rname for rname in unique(residues.names) if rname not in known_ligands]
        if len(rnames) == 0:
            self._log ( " - all residues have params - 1" )
            return

        log = residues[0].structure.session.logger

        # get ligand params from param folder
        from os import path
        dir_path = path.dirname(path.realpath(__file__))
        param_dir = path.join ( dir_path, 'param' )

        # Load GAFF atom types needed by small molecules parameterizations
        if not hasattr(self, '_gaff_types_added'):
            # Moriarty and Case parameterization uses GAFF atom types.
            from os import path
            gaff_types = path.join(param_dir, 'gaff2.xml')
            forcefield.loadFile(gaff_types)
            #self._gaff_types_added = True

        rnames_not_found = []
        for rname in rnames :
            pfile = path.join ( param_dir, '%s.xml' % rname.upper() )
            self._log ( " - looking for %s" % pfile )
            try :
                forcefield.loadFile(pfile)
            except :
                self._log ( " - not loaded/found" )
                rnames_not_found2.append(rname)
                continue
            log.info ( "Found params for %s" % rname )


        rnames = rnames_not_found
        if len(rnames) == 0:
            self._log ( " - all residues have params - 2" )
            return

        # check for more params in Tug Atoms plugin...
        param_dir = self._ligand_parameters_directory()
        if param_dir != None :
            log.info ( "Tug Ligands plugin not found - install for more params" )
        if param_dir != None :
            rnames_not_found = []
            from os import path
            mc_ligand_zip = path.join(param_dir, 'moriarty_and_case_ligands.zip')
            self._log ( " -> %s " % mc_ligand_zip )
            from zipfile import ZipFile
            with ZipFile(mc_ligand_zip) as zf:
                for rname in rnames:
                    try:
                        with zf.open(rname+'.xml') as xml_file:
                            forcefield.loadFile(xml_file)
                    except KeyError:
                        rnames_not_found.append(rname)
                        continue
                    log.info ( "Found params in Tug Ligands for %s" % rname )

        rnames = rnames_not_found
        if len(rnames) == 0:
            self._log ( " - all residues have params - 3" )
            return

        # Warn about ligands without parameters
        log.warning('Could not find OpenMM parameters for %d residues: %s'
                    % (len(rnames), ', '.join(rnames)))



    def _ligand_parameters_directory(self):
        try:
            from chimerax import tug_ligands
        except ImportError:
            return None
        return tug_ligands.parameters_directory


    def _system_from_prmtop(self, prmtop_path, impcrd_path):
        # load in Amber input files
        prmtop = app.AmberPrmtopFile(prmtop_path)
        inpcrd = app.AmberInpcrdFile(incrd_path)
        from numpy import array, float64
        from openmm import unit
        positions = 10*array(inpcrd.positions.value_in_unit(unit.nanometers), float64)  # Angstroms

        # Ordered atoms to match inpcrd order.
        # PDB can have atoms for a chain not contiguous, e.g. chain A hetatm at end of file.
        # But inpcrd has reordered so all chain atoms are contiguous.
        atom_pos = atoms.scene_coords
        from chimerax.geometry import find_closest_points
        i1,i2,near = find_closest_points(positions, atom_pos, 1.0)
        from numpy import empty, int32
        ai = empty((len(i1),), int32)
        ai[i1] = near
        self.atoms = oatoms = atoms[ai]
#        diff = oatoms.scene_coords - positions
#        print ('diff', diff.max(), diff.min())
#        p2a = {tuple(int(x) for x in a.scene_coord):a for a in atoms}
#        oatoms = [p2a[tuple(int(x) for x in p)] for p in positions]
#        print('\n'.join('%s %s' %(str(a1), str(a2)) for a1,a2 in zip(atoms, oatoms)))
#        print ('inpcrd', positions[:10])
#        print ('atoms', atom_pos[:10])

        # prepare system and integrator
        system = prmtop.createSystem(nonbondedMethod=app.CutoffNonPeriodic,
                                     nonbondedCutoff=self._nonbonded_cutoff*unit.angstrom,
                                     constraints=app.HBonds,
                                     rigidWater=True,
                                     ewaldErrorTolerance=0.0005
        )
        self._system = system

    def _create_system(self, forcefield):
        # We first try with ignoreExternalBonds = False.  This will fail any amino acid chain
        # ends do not have proper end-capping atoms, a common situation when there are missing
        # segments.  The addh command does not add termination atoms for missing segments.
        #
        # Then we try creating the system with ignoreExternalBonds = True.  This fails if
        # disulfides between cysteines are present with amber14 forcefield this giving a multiple
        # matching templates for CYS (CYM, CYX) error probably because the cysteine sulfur is
        # missing an attached atom.
        #
        # constraints = HBonds means the length of covalent bonds involving
        # hydrogen atoms are fixed to allow larger integration time steps.
        # Constraints are not supported to atoms with particle mass = 0 which
        # indicates a fixed atom position, and will generate errors.
        from openmm import app
        from openmm import unit
        try:
            system = forcefield.createSystem(self._topology,
                                             nonbondedMethod=app.CutoffNonPeriodic,
                                             nonbondedCutoff=self._nonbonded_cutoff*unit.angstrom,
                                             constraints=app.HBonds,
                                             rigidWater=True,
                                             ignoreExternalBonds=False)
        except Exception as e1:
            system = None
            err1 = e1
            self._mol.session.logger.warning ( "Forcefield error with EB: %s" % str(err1) )
            self._log ( "Forcefield error with EB: %s" % str(err1) )

        if system is not None:
            return system, None

        prob_res = None
        try:
            system = forcefield.createSystem(self._topology,
                                             nonbondedMethod=app.CutoffNonPeriodic,
                                             nonbondedCutoff=self._nonbonded_cutoff*unit.angstrom,
                                             constraints=app.HBonds,
                                             rigidWater=True,
                                             ignoreExternalBonds=True)
        except Exception as e2:
            #raise ForceFieldError('Missing atoms or parameterization needed by force field.\n' +
            #                      'All heavy atoms and hydrogens with standard names are required.\n' +
            #                      '\nError with ignoreExternalBonds=False was\n' + str(err1) +
            #                      '\nError with ignoreExternalBonds=True was\n' + str(e2))

            err2 = str(e2)
            self._log ( "Forcefield error with no EB: %s" % err2 )

            self._mol.session.logger.warning ( "Forcefield error no EB: %s" % err2 )
            if "Multiple non-identical matching templates found for residue" in err2 :
                ts = err2.split()
                resNum, resName = None, None
                for t in ts :
                    try :
                        resNum = int(t)
                    except :
                        pass
                    if t.startswith('(') and t.endswith('):') :
                        resName = t.replace('(','').replace('):','')
                if resNum != None and resName != None :
                    print ( " - problem res %d %s / %d res in %s" % (resNum, resName, len(self._mol.residues), self._mol.name) )
                    if resNum-1 >= 0 and resNum-1 < len( self._mol.residues ) :
                        res = self._mol.residues[resNum-1]
                        print ( " - found res number %d chain %s -- %s" % (res.number, res.chain_id, res.name) )
                        prob_res = res
            if "No template found for residue" in err2 :
                err2 = err2.split ( "." )[0]
                ts = err2.split()
                resNum, resName = None, None
                for t in ts :
                    try :
                        resNum = int(t)
                    except :
                        pass
                    if t.startswith('(') and t.endswith(')') :
                        resName = t.replace('(','').replace(')','')
                if resNum != None and resName != None :
                    print ( " - problem res %d %s / %d res in %s" % (resNum, resName, len(self._mol.residues), self._mol.name) )
                    if resNum-1 >= 0 and resNum-1 < len( self._mol.residues ) :
                        res = self._mol.residues[resNum-1]
                        print ( " - found res number %d chain %s -- %s" % (res.number, res.chain_id, res.name) )
                        prob_res = res


        return system, prob_res



    def getResNum ( rNum ) :
        for ri, res in enumerate ( self._mol.residues ) :
            if resNum != None :
                if resNum == ri+1 :
                    print ( "%s - %d - %d - %s" % (res.chain_id, ri+1, res.number, res.name) )
            else :
                print ( "%s - %d - %d - %s" % (res.chain_id, ri+1, res.number, res.name) )


    def _add_hydrogens(self, openmm_pdb, forcefield):
        # Need hydrogens to run simulation and most structures don't have them.
        # Most PDBs have partial residues (e.g. 1a0m with ARG with only CB) and adding hydrogens
        # give error of missing template in that case.  Eric says Chimera 1 uses rotamers tool
        # to extend the residues to include all heavy atoms.
        # When openmm does have all heavy atoms (1gcn) it preserves order of those atoms
        # but hydrogens are inserted after the heavy atom they are attached to.  Would need to
        # take that into account to map between ChimeraX structure and openmm structure with hydrogens.
        from openmm.app.modeller import Modeller
        m = Modeller(openmm_pdb.topology, openmm_pdb.positions)
        m.addHydrogens(forcefield)
        top, pos = m.getTopology(), m.getPositions()
        print('Before adding hydrogens')
        dump_topology(openmm_pdb.topology)
        print('After adding hydrogens')
        dump_topology(top)

class ForceFieldError(Exception):
    pass


def openmm_topology(atoms, bonds):
    '''Make OpenMM topology from ChimeraX atoms and bonds.'''
    a = atoms
    n = len(a)
    r = a.residues
    aname = a.names
    ename = a.element_names
    rname = r.names
    rnum = r.numbers
    cids = r.chain_ids
    from openmm.app import Topology, Element
    top = Topology()
    cmap = {}
    rmap = {}
    atoms = {}
    for i in range(n):
        cid = cids[i]
        if not cid in cmap:
            cmap[cid] = top.addChain()	# OpenMM chains have no name.
        rid = (rname[i], rnum[i], cid)
        if not rid in rmap:
            rmap[rid] = top.addResidue(rname[i], cmap[cid])
        element = Element.getBySymbol(ename[i])
        atoms[i] = top.addAtom(aname[i], element, rmap[rid])
    a1, a2 = bonds.atoms
    for i1, i2 in zip(a.indices(a1), a.indices(a2)):
        top.addBond(atoms[i1], atoms[i2])
    return top


_openmm_initialized = False
def initialize_openmm():
    # On linux need to set environment variable to find plugins.
    # Without this it gives an error saying there is no "CPU" platform.
    global _openmm_initialized
    if not _openmm_initialized:
        _openmm_initialized = True
        from sys import platform
        if platform == 'linux' or platform == 'darwin':
            from os import environ, path
            from chimerax import app_lib_dir
            environ['OPENMM_PLUGIN_DIR'] = path.join(app_lib_dir, 'plugins')


def dump_topology(t):
    for a in t.atoms():
        an, r, ai = a.name, a.residue, a.index
        rn, ri, c = r.name, r.index, r.chain
        ci = c.index
        print ('%d %s %s %d %d %s' % (ai, an, rn, ri, ci, c.id))

















# that's it

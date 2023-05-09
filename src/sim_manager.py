

##
##@Author: Greg Pintilie
##@Date:   14-Apr-2023
##@Email:  gregdp@gmail.com
##@Last modified by:   gregdp
##@Last modified time: 14-Apr-2023
##@License: Free for non-commercial use (see license.pdf)
##@Copyright: Copyright 2023 Grigore Pintilie
##





class SimManager () :

    def __init__ (self, session) :

        self._framesCaller = None
        self.session = session
        self._sim = None

        session.logger.info ( "Initialized sim manager" )



    def SetTemp ( self, temp ) :
        if self._sim != None :
            self._sim.SetTemp ( temp )


    def Minimize ( self, atoms, dmap, force_constant=10000, cutoff=10,
                    temperature=50, error_tolerance = 0.001,
                    steps=50, frames=50, finish=False ) :

        models = atoms.unique_structures
        print ( " - %d models " % len(models) )
        if len(models) > 1:
            from chimerax.core.errors import UserError
            raise UserError('More than one models (%d) specified' % len(us))

        mol = models[0]

        if dmap != None :
            print ( " - sim map set: %s" % dmap.name )

        from .sim import MolSim, ForceFieldError
        try:
            self._sim = MolSim ( atoms, dmap=dmap, force_constant = force_constant,
                                     cutoff = cutoff, temperature = temperature,
                                     tolerance = error_tolerance, steps = steps)
        except ForceFieldError as e:
            # Model could not be parameterized.
            from chimerax.core.errors import UserError
            raise UserError(str(e))

        self._sim.minimize ()


    def StartSim ( self, atoms, dmap=None, force_constant=10000, cutoff=10,
                    temperature=50, error_tolerance = 0.001,
                    steps=50, frames=50, finish=False ) :

        print ( " -- starting sim -- " )
        if self._sim :
            session.logger.error ( "Already running a simulation, stop first..." )
            return

        models = atoms.unique_structures
        print ( " - %d models " % len(models) )
        if len(models) > 1:
            from chimerax.core.errors import UserError
            raise UserError('More than one models (%d) specified' % len(us))

        mol = models[0]

        if dmap != None :
            print ( " - sim map set: %s" % dmap.name )


        from .sim import MolSim, ForceFieldError
        try:
            self._sim = MolSim ( atoms, dmap=dmap, force_constant = force_constant,
                                     cutoff = cutoff, temperature = temperature,
                                     tolerance = error_tolerance, steps = steps)
        except ForceFieldError as e:
            # Model could not be parameterized.
            from chimerax.core.errors import UserError
            raise UserError(str(e))

        #simo.tug_atoms ( atoms )
        #points = to_atoms.scene_coords
        self._sim.start ()
        self.session._davinci_sim = self._sim

        def simulation_frame ( session, frame ):
            #tugger.tug_to_positions(points)
            if self._sim :
                self._sim.step ()

        from chimerax.core.commands.motion import CallForNFrames
        numFrames = CallForNFrames.Infinite # 5000
        self._framesCaller = CallForNFrames ( simulation_frame, numFrames, self.session )

        if finish:
            from chimerax.std_commands.wait import wait
            wait(self.session, frames)


        #mm.add_mode(DavAtomsMode(session))
        print ( " - setting mousemode!" )
        from .sim_mm_pull import PullAtomMode
        mode = PullAtomMode ( self.session )
        mode._sim = self._sim
        self._prev_right_mode = self.session.ui.mouse_modes.mode ( 'right', [] )
        self.session.ui.mouse_modes.bind_mouse_mode('right', [], mode)



    def StopSim ( self ) :

        if self._framesCaller != None :
            print ( " - stopping after %d frames" % self._framesCaller.frame )
            self._framesCaller.done()
            self._framesCaller = None
        else :
            self.session.logger.info ( " - simulation not running" )

        if hasattr ( self, '_prev_right_mode' ) :
            self.session.ui.mouse_modes.bind_mouse_mode('right', [], self._prev_right_mode)
            print ( " - set previous right mouse mode" )













#

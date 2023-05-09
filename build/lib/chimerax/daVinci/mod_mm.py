

##
##@Author: Greg Pintilie
##@Date:   14-Apr-2023
##@Email:  gregdp@gmail.com
##@Last modified by:   gregdp
##@Last modified time: 14-Apr-2023
##@License: Free for non-commercial use (see license.pdf)
##@Copyright: Copyright 2023 Grigore Pintilie
##


write_logs = True

def register_mousemode(session):
    mm = session.ui.mouse_modes
    mm.add_mode ( PullAtomMode(session) )


from chimerax.mouse_modes import MouseMode
from .sim_mm_pull import Logger

class ModMouseMode ( MouseMode ) :

    name = 'davinci_model'
    icon_file = 'mod.png'

    def __init__(self, session):

        MouseMode.__init__(self, session)

        self._active = False
        self._handler = None
        self._last_frame_number = None
        self._log = Logger('/Users/greg/Desktop/mod.log' if write_logs else None)


    def mouse_down(self, event):
        #self._log('In mouse_down')
        x,y = event.position()
        self._log ( " - mouse down at %d %d" % (x,y) )
        print ( " - mouse down at %d %d" % (x,y) )
        view = self.session.main_view
        pick = view.picked_object(x,y)

        #self._pick_atom(pick)
        self._pick_map(pick)

        self.session._dav_pick = pick


    def _pick_map ( self, pick ) :
        if hasattr(pick, 'map') :
            P = pick.position
            self.modMan.AddAtPos ( P )


    def _pick_atom(self, pick):
        if hasattr(pick, 'atom'):
            a = pick.atom
            if not self._sim :
                self._log ( " - trying to pull but no simulation running ..." )
                return
            self._atom = a
            atom_i = self._sim.atoms.index ( a )
            self._log ( " - picked atom %s - %d to pull" % (a.name, atom_i) )
            self._active = True
            self._sim.pull_atom_i = atom_i


    def mouse_drag(self, event):
        x,y = event.position()
        #self._log ( " - mouse drag at %d %d" % (x,y) )
        #self._puller = Puller2D(x,y)
        #self._continue ()


    def mouse_up(self, event = None):

        self._log ( " - mouse up..." )
        return

        self._active = False
        self._last_frame_number = None

        if self._handler:
            self._log ( " - removed handler" )
            self.session.triggers.remove_handler(self._handler)
            self._handler = None
            amod = self._arrow_model
            if amod and not amod.deleted:
                amod.display = False
                amod.delete()
                self._arrow_model = None

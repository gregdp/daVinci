

##
##@Author: Greg Pintilie
##@Date:   14-Apr-2023
##@Email:  gregdp@gmail.com
##@Last modified by:   gregdp
##@Last modified time: 14-Apr-2023
##@License: Free for non-commercial use (see license.pdf)
##@Copyright: Copyright 2023 Grigore Pintilie
##


def register_mousemode(session):
    mm = session.ui.mouse_modes
    mm.add_mode ( PullAtomMode(session) )


from chimerax.mouse_modes import MouseMode

class PullAtomMode (MouseMode) :

    name = 'davinci_pull'
    icon_file = 'pull.png'

    def __init__(self, session):

        MouseMode.__init__(self, session)

        self._active = False
        self._handler = None
        self._last_frame_number = None
        self._puller = None
        self._arrow_model = None
        self._sim = None

        from os.path import exists
        logf = Logger.logFile() if exists ( Logger.sigFile() ) else None
        self._log = Logger( logf, "mm" )


    def mouse_down(self, event):
        #self._log('In mouse_down')
        x,y = event.position()
        self._log ( " - mouse down at %d %d" % (x,y) )
        view = self.session.main_view
        pick = view.picked_object(x,y)
        self._pick_atom(pick)


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
        self._puller = Puller2D(x,y)
        self._continue ()


    def mouse_up(self, event = None):

        self._log ( " - mouse up..." )
        self._active = False
        self._atom = None
        self._last_frame_number = None
        self._sim.pull_atom_i = None
        self._sim.pull_atom = False
        self._sim.reset_forces()

        if self._handler:
            self._log ( " - removed handler" )
            self.session.triggers.remove_handler(self._handler)
            self._handler = None
            amod = self._arrow_model
            if amod and not amod.deleted:
                amod.display = False
                amod.delete()
                self._arrow_model = None



    def _doit (self) :

        v = self.session.main_view
        if v.frame_number == self._last_frame_number:
            return	# not a new frame yet...

        self._last_frame_number = v.frame_number

        a = self._atom
        atom_xyz, target_xyz = self._puller.pull_to_point(a)

        #target_xyz.reshape((1,3))
        #self._log ( " -- pulling atom atpt -- %.3f %.3f %.3f" % (atom_xyz[0], atom_xyz[1], atom_xyz[2]) )
        #self._log ( " -- pulling target pt -- %.3f %.3f %.3f" % (target_xyz[0], target_xyz[1], target_xyz[2]) )

        self._log ( " - new - " )
        self._sim.pull_atom_pt = p = target_xyz
        self._sim.pull_atom = True
        self._log ( " - set pull atom %d pt %f %f %f" % (self._sim.pull_atom_i, p[0], p[1], p[2]) )

        atom_xyz, target_xyz = self._puller.pull_to_point(a)
        self._draw_arrow(target_xyz, atom_xyz)


    def _continue (self, *_):

        self._log ( " - mm continue ..." )

        if not self._active :
            return

        self._doit ()

        if self._handler is None:
            self._log ( " - adding new frame handler ..." )
            self._handler = self.session.triggers.add_handler ('new frame', self._continue)



    def _draw_arrow(self, xyz1, xyz2, radius = 0.1):
        self._log ( " - draw arrow" )
        a = self._arrow_model
        if a is None or a.deleted:
            from chimerax.core.models import Model
            s = self.session
            self._arrow_model = a = Model('Pull arrow', s)
            from chimerax.surface import cone_geometry
            v,n,t = cone_geometry(points_up = False)
            a.set_geometry(v, n, t)
            a.color = (0,255,0,255)
            s.models.add([a])
        # Scale and rotate prototype cylinder.
        from chimerax.atomic import structure
        from numpy import array, float32
        p = structure._bond_cylinder_placements(array(xyz1).reshape((1,3)),
                                                array(xyz2).reshape((1,3)),
                                                array([radius],float32))
        a.position = p[0]
        a.display = True

    def vr_press(self, event):
        # Virtual reality hand controller button press.
        view = self.session.main_view
        pick = event.picked_object(view)
        self._pick_atom(pick)

    def vr_motion(self, event):
        # Virtual reality hand controller motion.
        self._puller = Puller3D(event.tip_position)
        self._continue ()

    def vr_release(self, release):
        # Virtual reality hand controller button release.
        self.mouse_up()


class Puller2D:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def pull_to_point(self, atom):
        v = atom.structure.session.main_view
        x0,x1 = v.clip_plane_points(self.x, self.y)
        axyz = atom.scene_coord
        # Project atom onto view ray to get displacement.
        dir = x1 - x0
        from chimerax.geometry import inner_product
        pxyz = x0 + (inner_product(axyz - x0, dir)/inner_product(dir,dir)) * dir
        return axyz, pxyz


class Puller3D:
    def __init__(self, xyz):
        self.xyz = xyz

    def pull_to_point(self, atom):
        axyz = atom.scene_coord
        return axyz, self.xyz




class Logger:

    def __init__(self, filename = None, label = ""):
        self.filename = filename
        self._log_file = None
        self._label = label
        self._log_counter = 0

    @staticmethod
    def logFile () :
        return '/Users/greg/Desktop/dav.log'

    @staticmethod
    def sigFile () :
        return "/Users/greg/Dropbox/_mol/x/daVinci/"

    def __call__(self, message):
        if self.filename is None:
            return	# No logging
        self._log_file = f = open(self.filename,'a')
        f.write ( '%s_%04d ' % (self._label, self._log_counter) )
        f.write(message)
        f.write("\n")
        f.flush()
        self._log_counter += 1
        f.close()

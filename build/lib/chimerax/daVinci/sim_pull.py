

##
##@Author: Grigore Pintilie
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
    mm.add_mode(PullAtomMode(session))


from chimerax.mouse_modes import MouseMode

class PullAtomMode (MouseMode) :

    name = 'davinci_pull'
    icon_file = 'tug_atom.png'

    def __init__(self, session):

        MouseMode.__init__(self, session)

        self._active = False
        self._handler = None
        self._atom = None
        self._atom_i = None
        self._last_frame_number = None
        self._puller = None
        self._arrow_model = None
        self._sim = None

        self._log = Logger('/Users/greg/Desktop/mm.log' if write_logs else None)


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
            self._atom_i = ai = self._sim.atoms.index ( a )
            self._log ( " - picked atom %s - %d to pull" % (a.name, ai) )
            self._active = True


    def mouse_drag(self, event):
        x,y = event.position()
        #self._log ( " - mouse drag at %d %d" % (x,y) )
        self._puller = Puller2D(x,y)
        self._continue ()


    def mouse_up(self, event = None):

        self._log ( " - mouse up..." )
        self._active = False
        self._last_frame_number = None

        if self._handler:
            self._log ( " - removed handler" )
            self.session.triggers.remove_handler(self._handler)
            self._handler = None
            amod = self._arrow_model
            if amod and not amod.deleted:
                amod.display = False

        if self._atom_i != None :
            if self._sim._force != None :
                f = self._sim._force
                self._log ( " - stop - pull force on atom %s - %d" % (self._atom.name, self._atom_i) )
                f.setParticleParameters(self._atom_i, self._atom_i, (0,0,0,0))
                f.updateParametersInContext(self._sim._simulation.context)
                self._atom_i = None
                self._atom = None



    def _doit (self) :

        v = self.session.main_view
        if v.frame_number == self._last_frame_number:
            return	# not a new frame yet...

        self._last_frame_number = v.frame_number

        a = self._atom
        atom_xyz, target_xyz = self._puller.pull_to_point(a)

        #target_xyz.reshape((1,3))
        self._log ( " -- pulling atom atpt -- %.3f %.3f %.3f" % (atom_xyz[0], atom_xyz[1], atom_xyz[2]) )
        self._log ( " -- pulling target pt -- %.3f %.3f %.3f" % (target_xyz[0], target_xyz[1], target_xyz[2]) )

        if 1 :
            self._sim.pull_atom_to_pt ( self._atom, self._atom_i, target_xyz )

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
    def __init__(self, filename = None):
        self.filename = filename
        self._log_file = None
    def __call__(self, message, close = False):
        if self.filename is None:
            return	# No logging
        f = self._log_file
        if f is None:
            self._log_file = f = open(self.filename,'w')
            self._log_counter = 0
        f.write('%04d ' % self._log_counter)
        f.write(message)
        f.write("\n")
        f.flush()
        self._log_counter += 1
        if close:
            f.close()
            self._log_file = None

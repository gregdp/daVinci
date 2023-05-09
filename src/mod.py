

##
##@Author: Greg Pintilie
##@Date:   14-Apr-2023
##@Email:  gregdp@gmail.com
##@Last modified by:   gregdp
##@Last modified time: 14-Apr-2023
##@License: Free for non-commercial use (see license.pdf)
##@Copyright: Copyright 2023 Grigore Pintilie
##




class ModManager () :


    def __init__ ( self, session ) :

        print ( " - initialized ModManager" )
        self.session = session


    def AddRes ( self, resName ) :

        import chimerax

        toMod = None
        for m in self.session.models :
            if m.visible and type(m) == chimerax.atomic.structure.AtomicStructure :
                toMod = m

        if toMod == None :
            self.session.logger.error ( "Make sure a model is visible to add to" )
            return

        print ( "Adding to %s" % toMod.name )
        self.toMod = toMod
        self.resName = resName

        #mm.add_mode(DavAtomsMode(session))
        from .mod_mm import ModMouseMode
        mode = ModMouseMode ( self.session )
        mode.modMan = self

        print ( " - setting mousemode for adding water" )
        self._prev_right_mode = self.session.ui.mouse_modes.mode ( 'right', [] )
        self.session.ui.mouse_modes.bind_mouse_mode('right', [], mode)



    def AddAtPos ( self, P ) :

        print ( " - adding %s at  %.3f %.3f %.3f to %s" % (self.resName, P[0], P[1], P[2], self.toMod.name) )

        from chimerax.atomic.struct_edit import add_atom
        from chimerax.atomic import Sequence

        mol = self.toMod

        if self.resName.lower() == "w" or self.resName.lower() == "hoh" :

            last_res_num = 0
            last_res = None
            for res in mol.residues :
                if res.number > last_res_num :
                    last_res_num = res.number
                    last_res = res

            if last_res != None :
                isProt = Sequence.protein3to1(last_res.name) != "X"
                isNA = Sequence.nucleic3to1(res.name) != "X"
                if isProt or isNA :
                    print ( " - last (prot/na) res # %d" % last_res_num )
                    last_res += 10.0

            if last_res.number == last_res_num :
                print ( " - last (other) res # %d" % last_res_num )
                last_res_num += 1

            print ( " - adding at # %d" % last_res_num )

            newRes = mol.new_residue( "HOH", "A", last_res_num )
            P = mol.position.inverse() * P
            aO = add_atom ( "O", "O", newRes, P )
            aO.draw_mode = 2
            aO.color = (255, 13, 13, 255)

            from random import random
            from chimerax.geometry import normalize_vector, cross_product, place
            from chimerax.geometry import matrix, vector
            from numpy import array

            v1 = array ( [0,0,0] )
            while vector.length ( v1 ) < 1e-4 :
                v1 = array( [random()-0.5, random()-0.5, random()-0.5] )
            v1 = normalize_vector ( v1 )
            v1 = v1 * 1.0 # typical distance of O to H is 1.0

            v2 = array ( [0,0,0] )
            while vector.length ( v2 ) < 0.1 or vector.length ( vector.cross_product(v2,v1) ) < 0.1 :
                v2 = array( [random()-0.5, random()-0.5, random()-0.5] )
            v2 = normalize_vector ( v2 )
            ax = vector.cross_product ( v2, v1 )
            ax = normalize_vector ( ax )

            #t0 = matrix.translation_matrix ( P * -1 )
            #r = matrix.rotation_from_axis_angle ( ax, random() * 104.5 )
            #v2 = matrix.apply_matrix ( r, v )

            R = place.rotation ( ax, 104.5, center=[0,0,0] )
            v2 = R * v1
            print ( v1 )
            print ( v2 )

            aH1 = add_atom ( "H", "H1", newRes, P + v1 )
            aH1.draw_mode = 2
            aH1.color = (255, 255, 255, 255)

            aH2 = add_atom ( "H", "H2", newRes, P + v2 )
            aH2.draw_mode = 2
            aH2.color = (255, 255, 255, 255)

            mol.new_bond(aO, aH1)
            mol.new_bond(aO, aH2)
            # structure.reorder_residues(list(prev_residues) + residues1 + list(reversed(residues2)))
















#

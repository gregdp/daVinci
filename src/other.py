

##
##@Author: Greg Pintilie
##@Date:   14-Apr-2023
##@Email:  gregdp@gmail.com
##@Last modified by:   gregdp
##@Last modified time: 14-Apr-2023
##@License: Free for non-commercial use (see license.pdf)
##@Copyright: Copyright 2023 Grigore Pintilie
##



import chimerax

def listRes ( session, resNum=None ) :

    for m in session.models :
        if m.visible == True and type ( m ) == chimerax.atomic.structure.AtomicStructure :
            print ( m.name )

            for ri, res in enumerate ( m.residues ) :
                if resNum != None :
                    if resNum == ri+1 :
                        print ( "%s - %d - %d - %s" % (res.chain_id, ri+1, res.number, res.name) )
                else :
                    print ( "%s - %d - %d - %s" % (res.chain_id, ri+1, res.number, res.name) )


def ModId ( mol ) :
    mid = "%d" % mol.id[0]
    for i in range ( 1, len(mol.id) ) :
        mid += ".%d" % mol.id[i]
    return mid

def showLigands ( session ) :

    from chimerax.atomic import Sequence

    for m in session.models :
        # print ( "-- %s --" % m.name, type(m) )
        if m.visible and type(m) == chimerax.atomic.structure.AtomicStructure :
            modId = ModId(m)
            for res in m.residues :
                isProt = Sequence.protein3to1(res.name) != "X"
                isNA = Sequence.nucleic3to1(res.name) != "X"
                if not isProt and not isNA :
                    ss = "#%s/%s:%s" % (modId, res.chain_id, res.number)
                    print ( "Chain %s:%d %s/%s -- %s" % (res.chain_id, res.number, res.name, Sequence.protein3to1(res.name), ss) )


def showSelH ( session ) :

    ress = session.selection.items ("residues")[0] # chimerax.atomic.molarray.Residues
    print ( " - %d selected residues" % len(ress) )



# dav res resNum 646
# swap /E:96 CYS criteria d density #3
# swap /E:96 CYS criteria c











#

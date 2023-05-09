
##
##@Author: Greg Pintilie
##@Date:   14-Apr-2023
##@Email:  gregdp@gmail.com
##@Last modified by:   gregdp
##@Last modified time: 14-Apr-2023
##@License: Free for non-commercial use (see license.pdf)
##@Copyright: Copyright 2023 Grigore Pintilie
##


from chimerax.core.commands import CmdDesc
import chimerax


def current_simm (session) :
    if hasattr ( session, 'daVinci_sim_manager' ) :
        return session.daVinci_sim_manager

def set_current_simm ( session, simm ) :
    session.daVinci_sim_manager = simm


def daVinciCmd ( session, op, atoms=None, temp=50.0, inMap=None, resNum=None, add="" ) :

        # All command functions are invoked with ``session`` as its
        # first argument.  Useful session attributes include:
        #   logger: chimerax.core.logger.Logger instance
        #   models: chimerax.core.models.Models instance

        log = session.logger.info

        log ( "cmd __%s__" % ( op ) )

        if op == "mod" :
            print ( "adding... %s" % add )
            from .mod import ModManager
            modman = ModManager ( session )
            modman.AddRes ( add )

        elif op == "sim" or op == "s" :
            DoSim ( session, atoms, temp, inMap )

        elif op == "minimize" or op == "m" or op == "min" :
            DoMin ( session, atoms, temp, inMap )

        elif op == "stop" :
            simm = current_simm ( session )
            if simm != None :
                simm.StopSim()

        elif op == "set" or op == "s" :
            print ( "cmd -- temp -- %.2f" % temp )
            simm = current_simm ( session )
            if simm != None :
                simm.SetTemp ( temp )

        elif op == "res" :
            #from .other import *
            listRes ( session, resNum=resNum )

        elif op == "showLigands" or op == "sl" :
            #atoms, dmap = GetAtomsAndMap ( session, atoms, inMap )
            #from .other import *
            showLigands ( session )

        else :
            session.logger.warning ( "cmd -- %s not known -- " % op )




def GetAtomsAndMap ( session, atoms, inMap ) :

    if atoms != None :
        if len(atoms) == 0 :
            raise DaVinciParamError('Could not start sim, no atoms found as given' )
            return
    else :
        session.logger.info ( "No atoms specified - using first visible model in session..." )
        for m in session.models :
            print ( " - %s -- " % m.name, type(m) )
            if m.visible and type(m) == chimerax.atomic.structure.AtomicStructure :
                atoms = m.atoms
                break
        if atoms == None or len(atoms) == 0 :
            #session.logger.error ( "Could not start sim, no model with atoms found" )
            raise DaVinciParamError('Could not start sim, no model with atoms found' )

    print ( " -- starting with %d atoms -- " % len(atoms) )

    dmap = None
    if inMap != None :
        print ( inMap, len(inMap) )
        for mod in session.models :
            if mod.id == inMap :
                print ( " - found map: %s" % mod.name )
                print ( type(mod) )
                if type(mod) == chimerax.map.volume.Volume :
                    dmap = mod
                else :
                    session.logger.error ( "inMap parameter does not seem to specify a map"  )
                    raise DaVinciParamError('inMap parameter does not seem to specify a map' )
                    return

        if dmap == None :
            sp = "" + inMap
            session.logger.error ( "Could not find a map as specified %s" % sp  )
            raise DaVinciParamError( "Could not find a map as specified %s" % sp )

        print ( " -- in map: %s" % dmap.name )

    simm = current_simm ( session )
    if simm != None :
        print ( "Stopping previous simulation" )
        simm.StopSim ()

    return atoms, dmap



class DaVinciParamError(Exception):
    pass


def DoMin ( session, atoms, temp, inMap ) :

    atoms, dmap = GetAtomsAndMap ( session, atoms, inMap )

    from .sim_manager import SimManager
    simm = SimManager ( session )
    #simm.StartSim ( atoms, dmap=dmap, temperature=temp )
    set_current_simm ( session, simm )
    simm.Minimize ( atoms, dmap )

    print ( "___________!" )


def DoSim ( session, atoms, temp, inMap ) :

    atoms, dmap = GetAtomsAndMap ( session, atoms, inMap )

    from .sim_manager import SimManager
    simm = SimManager ( session )
    simm.StartSim ( atoms, dmap=dmap, temperature=temp )
    set_current_simm ( session, simm )

    print ( "___________!" )



def register_davinci_command ( name, synopsis, logger ) :

    #tugcommand.register_tug_command(logger)

    if 0 :
        logger.info ( "register -- %s -- %s" % (name, synopsis) )
        desc = CmdDesc()
        desc.synopsis = synopsis
        from chimerax.core.commands import register
        register ( name, desc, daVinciCmd )

    if 1 :
        logger.info ( "register --**-- %s -- %s -- op --**-- " % (name, synopsis) )
        from chimerax.core.commands import register, StringArg, FloatArg, ModelIdArg, IntArg #, CenterArg, BoolArg
        from chimerax.atomic import AtomsArg
        desc = CmdDesc(required = [('op', StringArg)],
                       #keyword = [('to_atoms', AtomsArg), ('force_constant', FloatArg)],
                       keyword = [('temp', FloatArg),
                                    ('inMap', ModelIdArg),
                                    ('atoms', AtomsArg),
                                    ('resNum', IntArg),
                                    ('add', StringArg)],
                       required_arguments = ['op'],
                       synopsis=synopsis)
        register( name, desc, daVinciCmd, logger=logger)


if 0 :
    try :
        import sys
        session.logger.info ( " - arg:", sys.argv[0] )
        register_davinci_command ( "davinci", "davinci", session.logger )
        session.logger.info ( " - registered command" )
        print ( "done" )
    except :
        print ( "?" )
        pass


#register_davinci_command ()
#daVinci_desc = CmdDesc()


# toolshed uninstall davinci; devel build /Users/greg/Dropbox/_mol/x/daVinci/; devel install /Users/greg/Dropbox/_mol/x/daVinci/
#
# devel clean /Users/greg/Dropbox/_mol/x/daVinci/

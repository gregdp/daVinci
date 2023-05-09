
##
##@Author: Greg Pintilie
##@Date:   14-Apr-2023
##@Email:  gregdp@gmail.com
##@Last modified by:   gregdp
##@Last modified time: 14-Apr-2023
##@License: Free for non-commercial use (see license.pdf)
##@Copyright: Copyright 2023 Grigore Pintilie
##


from chimerax.core.toolshed import BundleAPI

class _DaVinciAPI(BundleAPI):

    api_version = 1

    @staticmethod
    def initialize(session, bundle_info):
        """Register steered md mouse mode."""
        print ( "__init__ initialize ***************************" )
        if session.ui.is_gui:
            from . import sim
            sim.register_mousemode(session)

    @staticmethod
    def finish(session, bundle_info):
        print ( "__init__ finish ***************************" )
        # TODO: remove mouse mode
        pass


    @staticmethod
    def register_command(bi, ci, logger) :

        # bi is an instance of chimerax.core.toolshed.BundleInfo
        # ci is an instance of chimerax.core.toolshed.CommandInfo
        # logger is an instance of chimerax.core.logger.Logger

        # This method is called once for each command listed
        # in bundle_info.xml.  Since we only listed one command,
        # we expect only a single call to this method.

        # We import the function to call and its argument
        # description from the ``cmd`` module, adding a
        # synopsis from bundle_info.xml if none is supplied
        # by the code.

        logger.info ( " - registering daVinci command - " )
        from . import cmd
        cmd.register_davinci_command(ci.name, ci.synopsis, logger)

        #desc = cmd.daVinci_desc
        #if desc.synopsis is None:
        #    desc.synopsis = ci.synopsis
        #    logger.info( "syn: %s" % desc.synopsis )

        #from chimerax.core.commands import register
        #register(ci.name, desc, cmd.daVinciCmd)

        if 0 :
            if logger.session.ui.is_gui:
                print ( " - registering mouse mode" )
                from . import sim
                sim.register_mousemode(logger.session)



bundle_api = _DaVinciAPI()

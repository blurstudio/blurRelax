import maya.cmds as cmds
import maya.mel as mel
import os

from PySide2 import QtGui

NAME_WIDGET = 'BlurRelax_name'
RADIUS_WIDGET = 'BlurRelax_radius'
NEW_BIND_MESH_WIDGET = 'BlurRelax_newbindmesh'
BIND_FILE_WIDGET = 'BlurRelax_bindfile'


BORDER_DDL = "BlurRelax_border_ddl"
HARD_DDL = "BlurRelax_hard_ddl"
GROUP_DDL = "BlurRelax_group_ddl"
REPROJ_CHK = "BlurRelax_reproj_chk"
REPROJ_SLD = "BlurRelax_reproj_sld"
PRESERVE_VOL_SLD = "BlurRelax_preserve_vol_sld"
ITERATIONS_SLD = "BlurRelax_iterations_sld"
DELTA_CHK = "BlurRelax_delta_chk"
DELTA_SLD = "BlurRelax_delta_sld"

WIDGET_NAMES = (
    BORDER_DDL,
    HARD_DDL,
    GROUP_DDL,
    REPROJ_CHK,
    REPROJ_SLD,
    PRESERVE_VOL_SLD,
    ITERATIONS_SLD,
    DELTA_CHK,
    DELTA_SLD,
)
MENU_ITEMS = []


def create_menuitems():
    global MENU_ITEMS
    if MENU_ITEMS:
        # Already created
        return
    if cmds.about(api=True) < 201600:
        cmds.warning('BlurRelax menus only available in Maya 2016 and higher.')
        return
    for menu in ['mainDeformMenu', 'mainRigDeformationsMenu']:
        # Make sure the menu widgets exist first.
        mel.eval('ChaDeformationsMenu MayaWindow|{0};'.format(menu))
        items = cmds.menu(menu, query=True, itemArray=True)
        for item in items:
            if cmds.menuItem(item, query=True, divider=True):
                section = cmds.menuItem(item, query=True, label=True)
            menu_label = cmds.menuItem(item, query=True, label=True)
            if menu_label == 'Delta Mush':
                if section == 'Create':
                    BlurRelax_item = cmds.menuItem(
                        label="BlurRelax",
                        command=create_blur_relax,
                        sourceType='python',
                        insertAfter=item,
                        parent=menu,
                    )
                    BlurRelax_options = cmds.menuItem(
                        command=display_blur_relax_options,
                        insertAfter=BlurRelax_item,
                        parent=menu,
                        optionBox=True,
                    )
                    MENU_ITEMS.append(BlurRelax_item)
                    MENU_ITEMS.append(BlurRelax_options)
            elif menu_label == 'Delta Mush' and section == 'Paint Weights':
                item = cmds.menuItem(
                    label="BlurRelax",
                    command=paint_blur_relax_weights,
                    sourceType='python',
                    insertAfter=item,
                    parent=menu,
                )
                MENU_ITEMS.append(item)


def create_blur_relax(*args, **kwargs):
    cmds.loadPlugin('BlurRelax', quiet=True)
    nodes = cmds.deformer(type="BlurRelax")
    kwargs = get_create_command_kwargs()
    for node in nodes:
        for attr, value in kwargs.items():
            cmds.setAttr("{0}.{1}".format(node, attr), value)


def get_create_command_kwargs():
    """Gets the BlurRelax command arguments either from the option box widgets or the saved
    option vars.  If the widgets exist, their values will be saved to the option vars.
    @return A dictionary of the kwargs to the BlurRelax command."""
    kwargs = {}

    if cmds.optionMenu(BORDER_DDL, exists=True):
        # The index is 1-based
        val = cmds.optionMenu(BORDER_DDL, query=True, select=True) - 1
        kwargs['borderBehavior'] = val
        cmds.optionVar(intValue=(BORDER_DDL, val))
    else:
        kwargs['borderBehavior'] = cmds.optionVar(query=BORDER_DDL)

    if cmds.optionMenu(HARD_DDL, exists=True):
        # The index is 1-based
        val = cmds.optionMenu(HARD_DDL, query=True, select=True) - 1
        kwargs['hardEdgeBehavior'] = val
        cmds.optionVar(intValue=(HARD_DDL, val))
    else:
        kwargs['hardEdgeBehavior'] = cmds.optionVar(query=HARD_DDL)

    if cmds.optionMenu(GROUP_DDL, exists=True):
        # The index is 1-based
        val = cmds.optionMenu(GROUP_DDL, query=True, select=True) - 1
        kwargs['groupEdgeBehavior'] = val
        cmds.optionVar(intValue=(GROUP_DDL, val))
    else:
        kwargs['groupEdgeBehavior'] = cmds.optionVar(query=GROUP_DDL)

    if cmds.checkBoxGrp(REPROJ_CHK, exists=True):
        val = cmds.checkBoxGrp(REPROJ_CHK, query=True, value1=True)
        kwargs['reproject'] = bool(val)
        cmds.optionVar(intValue=(REPROJ_CHK, int(val)))
    else:
        kwargs['reproject'] = bool(cmds.optionVar(query=REPROJ_CHK))

    if cmds.checkBoxGrp(DELTA_CHK, exists=True):
        val = cmds.checkBoxGrp(DELTA_CHK, query=True, value1=True)
        kwargs['delta'] = bool(val)
        cmds.optionVar(intValue=(DELTA_CHK, int(val)))
    else:
        kwargs['delta'] = bool(cmds.optionVar(query=DELTA_CHK))

    if cmds.floatSliderGrp(REPROJ_SLD, exists=True):
        val = cmds.floatSliderGrp(REPROJ_SLD, query=True, value=True)
        kwargs['reprojectDivs'] = val
        cmds.optionVar(floatValue=(REPROJ_SLD, val))
    else:
        kwargs['reprojectDivs'] = cmds.optionVar(query=REPROJ_SLD)

    if cmds.floatSliderGrp(PRESERVE_VOL_SLD, exists=True):
        val = cmds.floatSliderGrp(PRESERVE_VOL_SLD, query=True, value=True)
        kwargs['preserveVolume'] = val
        cmds.optionVar(floatValue=(PRESERVE_VOL_SLD, val))
    else:
        kwargs['preserveVolume'] = cmds.optionVar(query=PRESERVE_VOL_SLD)

    if cmds.floatSliderGrp(ITERATIONS_SLD, exists=True):
        val = cmds.floatSliderGrp(ITERATIONS_SLD, query=True, value=True)
        kwargs['iterations'] = val
        cmds.optionVar(floatValue=(ITERATIONS_SLD, val))
    else:
        kwargs['iterations'] = cmds.optionVar(query=ITERATIONS_SLD)

    if cmds.floatSliderGrp(DELTA_SLD, exists=True):
        val = cmds.floatSliderGrp(DELTA_SLD, query=True, value=True)
        kwargs['deltaMultiplier'] = val
        cmds.optionVar(floatValue=(DELTA_SLD, val))
    else:
        kwargs['deltaMultiplier'] = cmds.optionVar(query=DELTA_SLD)

    return kwargs



def display_blur_relax_options(*args, **kwargs):
    cmds.loadPlugin('BlurRelax', qt=True)
    layout = mel.eval('getOptionBox')
    cmds.setParent(layout)
    cmds.columnLayout(adj=True)

    for widget in WIDGET_NAMES:
        # Delete the widgets so we don't create multiple controls with the same name
        try:
            cmds.deleteUI(widget, control=True)
        except RuntimeError:
            pass

    if not cmds.optionVar(exists=BORDER_DDL):
        init_option_vars()

    pinBehaviors = ["None", "Pin", "Slide"]

    cmds.optionMenu(BORDER_DDL, label="Border Behavior")
    for b in pinBehaviors:
        cmds.menuItem(label=b, parent=BORDER_DDL)
    cmds.optionMenu(
        BORDER_DDL, edit=True, value=pinBehaviors[cmds.optionVar(query=BORDER_DDL)]
    )

    cmds.optionMenu(HARD_DDL, label="Hard Edge Behavior")
    for b in pinBehaviors:
        cmds.menuItem(label=b, parent=HARD_DDL)
    cmds.optionMenu(
        HARD_DDL, edit=True, value=pinBehaviors[cmds.optionVar(query=HARD_DDL)]
    )

    cmds.optionMenu(GROUP_DDL, label="Group Edge Behavior")
    for b in pinBehaviors:
        cmds.menuItem(label=b, parent=GROUP_DDL)
    cmds.optionMenu(
        GROUP_DDL, edit=True, value=pinBehaviors[cmds.optionVar(query=GROUP_DDL)]
    )

    cmds.checkBoxGrp(
        REPROJ_CHK,
        numberOfCheckBoxes=1,
        label='reproject',
        value1=cmds.optionVar(query=REPROJ_CHK),
    )

    cmds.floatSliderGrp(
        REPROJ_SLD,
        label='Reproject Divs',
        field=True,
        minValue=0,
        maxValue=3,
        fieldMinValue=0,
        fieldMaxValue=3,
        step=1,
        precision=0,
        value=cmds.optionVar(query=REPROJ_SLD),
    )

    cmds.floatSliderGrp(
        PRESERVE_VOL_SLD,
        label='Preserve Volume',
        field=True,
        minValue=0.0,
        maxValue=2.0,
        fieldMinValue=0.0,
        fieldMaxValue=2.0,
        step=0.01,
        precision=2,
        value=cmds.optionVar(query=PRESERVE_VOL_SLD),
    )

    cmds.floatSliderGrp(
        ITERATIONS_SLD,
        label='Iterations',
        field=True,
        minValue=0.0,
        maxValue=50.0,
        fieldMinValue=0.0,
        fieldMaxValue=1000.0,
        step=0.1,
        precision=1,
        value=cmds.optionVar(query=ITERATIONS_SLD),
    )

    cmds.checkBoxGrp(
        DELTA_CHK,
        numberOfCheckBoxes=1,
        label='Delta',
        value1=cmds.optionVar(query=DELTA_CHK),
    )
    cmds.floatSliderGrp(
        DELTA_SLD,
        label='Delta Multiplier',
        field=True,
        minValue=0.0,
        maxValue=1.0,
        fieldMinValue=0.0,
        fieldMaxValue=10.0,
        step=0.1,
        precision=1,
        value=cmds.optionVar(query=DELTA_SLD),
    )

    mel.eval('setOptionBoxTitle("BlurRelax Options");')
    mel.eval('setOptionBoxCommandName("BlurRelax");')
    apply_close_button = mel.eval('getOptionBoxApplyAndCloseBtn;')
    cmds.button(apply_close_button, edit=True, command=apply_and_close)
    apply_button = mel.eval('getOptionBoxApplyBtn;')
    cmds.button(apply_button, edit=True, command=create_blur_relax)
    reset_button = mel.eval('getOptionBoxResetBtn;')
    # For some reason, the buttons in the menu only accept MEL.
    cmds.button(
        reset_button,
        edit=True,
        command='python("import BlurRelax.menu; BlurRelax.menu.reset_to_defaults()");',
    )
    close_button = mel.eval('getOptionBoxCloseBtn;')
    cmds.button(close_button, edit=True, command=close_option_box)
    save_button = mel.eval('getOptionBoxSaveBtn;')
    cmds.button(
        save_button,
        edit=True,
        command='python("import BlurRelax.menu; BlurRelax.menu.get_create_command_kwargs()");',
    )
    mel.eval('showOptionBox')


def apply_and_close(*args, **kwargs):
    """Create the BlurRelax deformer and close the option box."""
    create_blur_relax()
    mel.eval('saveOptionBoxSize')
    close_option_box()


def close_option_box(*args, **kwargs):
    mel.eval('hideOptionBox')


def reset_to_defaults(*args, **kwargs):
    """Reset the BlurRelax option box widgets to their defaults."""
    cmds.optionMenu(BORDER_DDL, edit=True, value="Pin")
    cmds.optionMenu(HARD_DDL, edit=True, value="None")
    cmds.optionMenu(GROUP_DDL, edit=True, value="None")
    cmds.checkBoxGrp(REPROJ_CHK, edit=True, value1=False)
    cmds.checkBoxGrp(DELTA_CHK, edit=True, value1=False)
    cmds.floatSliderGrp(REPROJ_SLD, edit=True, value=1)
    cmds.floatSliderGrp(PRESERVE_VOL_SLD, edit=True, value=0.0)
    cmds.floatSliderGrp(ITERATIONS_SLD, edit=True, value=10.0)
    cmds.floatSliderGrp(DELTA_SLD, edit=True, value=1.0)


def init_option_vars():
    """Initialize the option vars the first time the ui is run"""
    cmds.optionVar(intValue=(BORDER_DDL, 1))
    cmds.optionVar(intValue=(HARD_DDL, 0))
    cmds.optionVar(intValue=(GROUP_DDL, 0))
    cmds.optionVar(intValue=(REPROJ_CHK, 0))
    cmds.optionVar(intValue=(DELTA_CHK, 0))
    cmds.optionVar(floatValue=(REPROJ_SLD, 1.0))
    cmds.optionVar(floatValue=(PRESERVE_VOL_SLD, 0.0))
    cmds.optionVar(floatValue=(ITERATIONS_SLD, 10.0))
    cmds.optionVar(floatValue=(DELTA_SLD, 1.0))


def get_wrap_node_from_object(obj):
    """Get a wrap node from the selected geometry."""
    if cmds.nodeType(obj) == 'BlurRelax':
        return obj
    history = cmds.listHistory(obj, pdo=0) or []
    wrap_nodes = [node for node in history if cmds.nodeType(node) == 'BlurRelax']
    if not wrap_nodes:
        raise RuntimeError('No BlurRelax node found on {0}.'.format(obj))
    if len(wrap_nodes) == 1:
        return wrap_nodes[0]
    else:
        # Multiple wrap nodes are deforming the mesh.  Let the user choose which one
        # to use.
        return QtGui.QInputDialog.getItem(
            None, 'Select BlurRelax node', 'blurRelax node:', wrap_nodes
        )


def get_wrap_node_from_selected():
    """Get a wrap node from the selected geometry."""
    sel = cmds.ls(sl=True) or []
    if not sel:
        raise RuntimeError('No BlurRelax found on selected.')
    return get_wrap_node_from_object(sel[0])


def destroy_menuitems():
    """Remove the BlurRelax items from the menus."""
    global MENU_ITEMS
    for item in MENU_ITEMS:
        cmds.deleteUI(item, menuItem=True)
    MENU_ITEMS = []


def paint_blur_relax_weights(*args, **kwargs):
    """Activates the paint BlurRelax weights context."""
    sel = cmds.ls(sl=True)
    if not sel:
        return
    wrap_node = get_wrap_node_from_selected()
    if not wrap_node:
        return
    mel.eval(
        'artSetToolAndSelectAttr("artAttrCtx", "BlurRelax.{0}.weights");'.format(
            wrap_node
        )
    )

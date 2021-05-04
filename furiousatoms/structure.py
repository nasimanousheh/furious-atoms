import vtk
import numpy as np
from fury import actor


def mobius():
    uResolution = 51
    vResolution = 51
    surface = vtk.vtkParametricMobius()
    surface.SetMinimumV(-0.25)
    surface.SetMaximumV(0.25)

    source = vtk.vtkParametricFunctionSource()
    source.SetUResolution(uResolution)
    source.SetVResolution(vResolution)
    source.SetParametricFunction(surface)
    source.Update()

    return source


def bbox(box_lx, box_ly, box_lz, colors=(0, 0, 0),
         linewidth=1, fake_tube=True, bbox_type='line'):
    """Return Bounding Box actor.

    Parameters
    ----------
    box_lx : int
        [description]
    box_ly : int
        [description]
    box_lz : int
        [description]
    colors : tuple, optional
        [description], by default (0, 0, 0)
    linewidth : int, optional
        [description], by default 1
    fake_tube : bool, optional
        [description], by default True
    bbox_type : str, optional
        [description], by default 'line'

    Returns
    -------
    [type]
        [description]
    """

    edge1 = 0.5 * np.array([[box_lx, box_ly, box_lz],
                            [box_lx, box_ly, -box_lz],
                            [-box_lx, box_ly, -box_lz],
                            [-box_lx, box_ly, box_lz],
                            [box_lx, box_ly, box_lz]])
    edge2 = 0.5 * np.array([[box_lx, box_ly, box_lz],
                            [box_lx, -box_ly, box_lz]])
    edge3 = 0.5 * np.array([[box_lx, box_ly, -box_lz],
                            [box_lx, -box_ly, -box_lz]])
    edge4 = 0.5 * np.array([[-box_lx, box_ly, -box_lz],
                            [-box_lx, -box_ly, -box_lz]])
    edge5 = 0.5 * np.array([[-box_lx, box_ly, box_lz],
                            [-box_lx, -box_ly, box_lz]])
    lines = [edge1, -edge1, edge2, edge3, edge4, edge5]

    if bbox_type == 'cube':
        box_centers = np.array([[0, 0, 25]])
        box_directions = np.array([[0, 1, 0]])
        bbox_actor = actor.box(box_centers, box_directions, colors,
                               scales=(box_lx, box_ly, box_lz))
        bbox_actor.GetProperty().SetRepresentationToWireframe()
        bbox_actor.GetProperty().SetLineWidth(linewidth)

    else:
        bbox_actor = actor.line(lines, colors=colors,
                                linewidth=linewidth,
                                fake_tube=fake_tube)

    return bbox_actor, lines


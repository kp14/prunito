import itertools
from math import acos, atan2, cos, exp, pi, sin, sqrt
from numpy import array, float64, zeros
from numpy.linalg import norm


def scalar_vector_product(first, second):
    """Compute product of cartesian vector coordinates.

    Args:
        self, first, second

    Returns:
        coordinates [x, y, z]
    """

    product = []
    for a, b in itertools.izip(first, second):
        product.append(a + b)
    return product


def logistic_sigmoid(x, a):
    """Computes so-called logistic curve, a sigmoidal function used in modelling
    population growth etc. In this implementation, the parameter a regulates the slope
    or 'growth rate' of the sigmoid's rising portion. When a=0, the function collapses
    to the identity function y=x.

    See: http://www.flong.com/texts/code/shapers_exp/

    Args:
        x, a

    Returns:
        y
    """
    x = float(x)
    epsilon = 0.0001
    min_param_a = 0.0 + epsilon
    max_param_a = 1.0 - epsilon
    a = max(min_param_a, min(max_param_a, a))
    a = (1 / (1 - a) - 1)

    A = 1.0 / (1.0 + exp(0 - ((x - 0.5) * a * 2.0)))
    B = 1.0 / (1.0 + exp(a))
    C = 1.0 / (1.0 + exp(0 - a))
    y = (A - B) / (C - B)
    return y


def remap(val, max_, slope=0):
    """Map a range of input data, the domain, to an output range
    using the logistic sigmoidal function for mapping.

    Args:
        val, max_, slope

    Returns:
        range boundaries as tuple

    """
    out_min = 0.0
    out_max = 1.0
    in_min = 0.0
    in_max = max_

    if val > in_max:
        val = in_max

    x = float(val / in_max)

    return logistic_sigmoid(x, slope)


def cartesian_to_spherical(vector):
    """Convert the Cartesian vector [x, y, z] to spherical coordinates [r, theta, phi].

    The parameter r is the radial distance, theta is the polar angle, and phi is the azimuth.


    @param vector:  The Cartesian vector [x, y, z].
    @type vector:   numpy rank-1, 3D array
    @return:        The spherical coordinate vector [r, theta, phi].
    @rtype:         numpy rank-1, 3D array
    """

    # The radial distance.
    r = norm(vector)

    # Unit vector.
    unit = vector / r

    # The polar angle.
    theta = acos(unit[2])

    # The azimuth.
    phi = atan2(unit[1], unit[0])

    # Return the spherical coordinate vector.
    return array([r, theta, phi], float64)


def spherical_to_cartesian(spherical_vect):
    """Convert the spherical coordinate vector [r, theta, phi] to the Cartesian vector [x, y, z].

    The parameter r is the radial distance, theta is the polar angle, and phi is the azimuth.


    @param spherical_vect:  The spherical coordinate vector [r, theta, phi].
    @type spherical_vect:   3D array or list

    """

    # Trig alias.
    sin_theta = sin(spherical_vect[1])

    # The vector.
    cart_vect = []
    cart_vect.append(spherical_vect[0] * cos(spherical_vect[2]) * sin_theta)  #x
    cart_vect.append(spherical_vect[0] * sin(spherical_vect[2]) * sin_theta)  #y
    cart_vect.append(spherical_vect[0] * cos(spherical_vect[1]))  #z

    return cart_vect


def get_points_equiangularly_distanced_on_sphere(numberOfPoints=10):
    """ each point you get will be of form 'x, y, z'; in cartesian coordinates
        eg. the 'l2 distance' from the origin [0., 0., 0.] for each point will be 1.0
        ------------
        converted from:  http://web.archive.org/web/20120421191837/http://www.cgafaq.info/wiki/Evenly_distributed_points_on_sphere )
    """
    dlong = pi * (3.0 - sqrt(5.0))  # ~2.39996323
    dz = 2.0 / numberOfPoints
    long = 0.0
    z = 1.0 - dz / 2.0
    ptsOnSphere = []
    for k in range(0, numberOfPoints):
        r = sqrt(1.0 - z * z)
        ptNew = (cos(long) * r, sin(long) * r, z)
        ptsOnSphere.append(ptNew)
        z = z - dz
        long = long + dlong
    return ptsOnSphere


def get_points_equiangularly_distanced_on_circle(number_of_points=0):
    """Cartesian coordinates for points on circle with roughly equal distances between them.
    Radius of circle is set to 1.

    Args:
        number_of_points=0

    Returns:
        list of coordniate tuples
    """
    pts_on_circle = []
    degrees = 359
    radius = 1
    if number_of_points == 0:
        return pts_on_circle
    else:
        angle = float(degrees / number_of_points)

    # if __name__ == '__main__':


#     ptsOnSphere = GetPointsEquiAngularlyDistancedOnSphere( 80)
#
#     #toggle True/False to print them
#     if( True ):
#         for pt in ptsOnSphere:  print( pt)
#
#     #toggle True/False to plot them
#     if(True):
#         from numpy import *
#         import pylab as p
#         import mpl_toolkits.mplot3d.axes3d as p3
#
#         fig=p.figure()
#         ax = p3.Axes3D(fig)
#
#         x_s=[];y_s=[]; z_s=[]
#
#         for pt in ptsOnSphere:
#             x_s.append( pt[0]); y_s.append( pt[1]); z_s.append( pt[2])
#
#         ax.scatter3D( array( x_s), array( y_s), array( z_s) )
#         ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
#         p.show()
#         #end

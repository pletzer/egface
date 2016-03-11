import scipy.integrate
import scipy.linalg
import numpy
import math
import re
from math import sin, cos, tan, asin, acos, atan, atan2, pi, exp, log, log10, e

class FormSetter:

    def __init__(self, order):
        """
        Constructor
        @param order order of the form (0, 1, 2, or 3)
        """
        # order of the form (0 = point, 1 = edge, 2 = face, 3 = cell)
        self.order = order

        # dictionary basis element => function, eg
        # {(0, 2): 'x**2'} represents x**2 * dx ^ dz
        self.terms = {}

        # basis element pattern, eg 'd\s*([xyz])\s*\^\s*d\s*([xyz])\s*'
        self.formPat = ''
        dxyzPat = '\s*d\s*([xyz])\s*'
        wedgePat = '\^'
        if order >= 1:
            self.formPat += dxyzPat
        if order >= 2:
            self.formPat += wedgePat + dxyzPat
        if order == 3:
            self.formPat += wedgePat + dxyzPat
        if order > 0:
            self.formPat = '\s*\*' + self.formPat

        # edge lengths
        self.dVerts = []

        # base vertex
        self.baseVert = None   

    def setExpression(self, expr):
        """
        Set form expression
        @param expr expression, eg x**2 * dx ^ dz + sin(pi*y)* dy ^ dx
        """
        if self.formPat == '':
            self.terms[''] = expr
            return

        coordName2Index = {'x': 0, 'y': 1, 'z': 2}
        # break the expresssion into independent terms
        continueFlag = True
        e = expr[:] # copy
        coords = []
        while continueFlag:
            m = re.search(self.formPat, e)
            if m:
                coords.append(tuple([coordName2Index[m.group(1 + o)] for o in range(self.order)]))
                e = re.sub(re.escape(m.group(0)), ',', e, 1) # replace the first occurrence only
                e = re.sub(r'\,$', '', e)
            else:
                continueFlag = False
        self.terms = dict(zip(coords, e.split(',')))

    def setElement(self, *vertices):
        """
        Set the cell element vertices
        @param vertices (ordert + 1)-tuple of vertices of the point, edge, face, or cell
        """
        assert len(vertices) == self.order + 1
        self.dVerts = [numpy.array(v) - numpy.array(vertices[0]) for v in vertices[1:]]
        self.baseVert = numpy.array(vertices[0])

    def funcXsi(self, xsi):
        """
        1-d integrand
        """
        return eval(self.funcXsiStr)

    def funcXsiEta(self, xsi, eta):
        """
        2-d integrand
        """
        return eval(self.funcXsiEtaStr)

    def funcXsiEtaZet(self, xsi, eta, zet):
        """
        3-d integrand
        """
        return eval(self.funcXsiEtaZetStr)

    def evaluate(self):
        """
        Evaluate (integrate) the form over the element
        @return value
        """
        res = 0

        if self.order == 0:
            # simple evaluation
            for k, v in self.terms.items():
                v = re.sub('x', 'self.baseVert[0]', v)
                v = re.sub('y', 'self.baseVert[1]', v)
                v = re.sub('z', 'self.baseVert[2]', v)
                res += eval(v)

        elif self.order == 1:
            # integration along a line
            for k, v in self.terms.items():
                k0, = k
                det = self.dVerts[0][k[0]]
                v = re.sub('x', '(self.baseVert[0] + xsi*self.dVerts[0][0])', v)
                v = re.sub('y', '(self.baseVert[1] + xsi*self.dVerts[0][1])', v)
                v = re.sub('z', '(self.baseVert[2] + xsi*self.dVerts[0][2])', v)
                self.funcXsiStr = v
                res += scipy.integrate.quad(self.funcXsi, 0., 1.)[0] * det

        elif self.order == 2:
            # integration along a triangle
            for k, v in self.terms.items():
                k0, k1 = k
                det = self.dVerts[0][k0] * self.dVerts[1][k1] - self.dVerts[0][k1] * self.dVerts[1][k0]
                v = re.sub('x', '(self.baseVert[0] + xsi*self.dVerts[0][0] + eta*self.dVerts[1][0])', v)
                v = re.sub('y', '(self.baseVert[1] + xsi*self.dVerts[0][1] + eta*self.dVerts[1][1])', v)
                v = re.sub('z', '(self.baseVert[2] + xsi*self.dVerts[0][2] + eta*self.dVerts[1][2])', v)
                self.funcXsiEtaStr = v
                def loEta(xsi):
                    return 0.0
                def hiEta(xsi):
                    return 1.0 - xsi
                res += scipy.integrate.dblquad(self.funcXsiEta, 0., 1., loEta, hiEta)[0] * det

        elif self.order == 3:
            # integration over tetrahedron
            for k, v in self.terms.items():
                k0, k1, k2 = k
                jac = [self.dVerts[0], self.dVerts[1], self.dVerts[2]]
                det = scipy.linalg.det(jac)
                v = re.sub('x', '(self.baseVert[0] + xsi*self.dVerts[0][0] + eta*self.dVerts[1][0] + zet*self.dVerts[2][0])', v)
                v = re.sub('y', '(self.baseVert[1] + xsi*self.dVerts[0][1] + eta*self.dVerts[1][1] + zet*self.dVerts[2][1])', v)
                v = re.sub('z[^et]', '(self.baseVert[2] + xsi*self.dVerts[0][2] + eta*self.dVerts[1][2] + zet*self.dVerts[2][2])', v)
                self.funcXsiEtaZetStr = v
                def loEta(xsi):
                    return 0.0
                def hiEta(xsi):
                    return 1.0 - xsi
                def loZet(xsi, eta):
                    return 0.0
                def hiZet(xsi, eta):
                    return 1.0 - xsi - eta
                res += scipy.integrate.tplquad(self.funcXsiEtaZet, 0., 1., loEta, hiEta, loZet, hiZet)[0] * det

        else:
            raise RunTimeError, 'ERROR: 0 <= order <= 3 but got {0}'.format(order)

        return res

###############################################################################
def test0():
    fs = FormSetter(0)
    fs.setExpression('x*y + z')
    fs.setElement([1., 2., 3.])
    res = fs.evaluate()
    assert math.fabs(res - (1.*2. + 3.)) < 1.e-10

def test1():
    fs = FormSetter(1)
    fs.setExpression('x* dy + z**2 *dz')
    fs.setElement([1., 2., 3.], [1.1, 2.2, 3.3])
    res = fs.evaluate()
    assert math.fabs(res - 0.21 - 2.979) < 1.e-10

def test2():
    fs = FormSetter(2)
    fs.setExpression('x* dy ^ dz')
    fs.setElement([1., 2., 3.], [1.1, 2.2, 3.3], [-0.1, 0.2, 0.5])
    res = fs.evaluate()
    assert math.fabs(res - 0.0133333333333) < 1.e-8

def test3():
    fs = FormSetter(3)
    fs.setExpression('x* dy ^ dz ^ dx')
    fs.setElement([1., 2., 3.], [1.1, 2.2, 3.3], [-0.1, 0.2, 0.5], [0.2, 0.3, 1.4])
    res = fs.evaluate()
    assert math.fabs(res - 0.00366666666667) < 1.e-8

def testDblQuad():
    xa, xb, xc = 1., 1.1, -0.1
    ya, yb, yc = 2., 2.2, 0.2
    za, zb, zc = 3., 3.3, 0.5
    dxb, dyb, dzb = xb - xa, yb - ya, zb - za
    dxc, dyc, dzc = xc - xa, yc - ya, zc - za
    def f(xsi, eta):
        return xa + xsi*(xb - xa) + eta*(xc - xa)
    def etaMin(xsi):
        return 0.0
    def etaMax(xsi):
        return 1.0 - xsi

    area = dyb * dzc - dyc * dzb

    res = scipy.integrate.dblquad(f, 0., 1., etaMin, etaMax)[0] * area
    print 'testDblQuad res = ', res

def testTplQuad():
    xa, xb, xc, xd = 1., 1.1, -0.1, 0.2
    ya, yb, yc, yd = 2., 2.2, 0.2, 0.3
    za, zb, zc, zd = 3., 3.3, 0.5, 1.4

    dxb, dyb, dzb = xb - xa, yb - ya, zb - za
    dxc, dyc, dzc = xc - xa, yc - ya, zc - za
    dxd, dyd, dzd = xd - xa, yd - ya, zd - za

    def f(xsi, eta, zet):
        return xa + xsi*(xb - xa) + eta*(xc - xa) + zet*(xd - xa)
    def etaMin(xsi):
        return 0.0
    def etaMax(xsi):
        return 1.0 - xsi
    def zetMin(xsi, eta):
        return 0.0
    def zetMax(xsi, eta):
        return 1.0 - xsi - eta

    volume = dyb*(dzc*dxd - dzd*dxc) - dyc*(dzb*dxd - dzd*dxb) + dyd*(dzb*dxc - dzc*dxb)
    vol2 = scipy.linalg.det([[dyb, dzb, dxb], [dyc, dzc, dxc], [dyd, dzd, dxd]])
    assert math.fabs(volume - vol2) < 1.e-10

    res = scipy.integrate.tplquad(f, 0., 1., etaMin, etaMax, zetMin, zetMax)[0] * volume
    print 'test TplQuad res = ', res

if __name__ == '__main__': 
    testDblQuad()
    testTplQuad()
    test0()
    test1()
    test2()
    test3()
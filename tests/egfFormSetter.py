import scipy.integrate
import numpy
import math
import re

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

        # basis element pattern, eg 'd\s*([xyz])\s*\^d\s*([xyz])\s*'
        self.formPat = re.sub(r'^\^', '', reduce(lambda x, y: x + '^' + y, ['d\s*([xyz])\s*']*order, ''))
        if self.formPat != '':
            self.formPat = '\s*\\*\s*' + self.formPat

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

    def evaluate(self):
        """
        Evaluate (integrate) the form over the element
        @return value
        """
        from math import sin, cos, tan, asin, acos, atan, atan2, pi, exp, log, log10, e
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
                def f(xsi, eta):
                    return eval(v)
                def loEta(xsi):
                    return 0.0
                def hiEta(xsi):
                    return 1.0 - xsi
                res += scipy.integrate.dblquad(self.funcXsiEta, 0., 1., loEta, hiEta)[0] * det

        elif self.order == 3:
            # integration over tetrahedron
            for k, v in self.terms.items():
                k0, k1, k2 = k
                det = scipy.linalg.det(jac)
                v = re.sub('x', '(self.baseVert[0] + xsi*self.dVerts[0][0] + eta*self.dVerts[1][0] + zet*self.dVerts[2][0])', v)
                v = re.sub('y', '(self.baseVert[1] + xsi*self.dVerts[0][1] + eta*self.dVerts[1][1] + zet*self.dVerts[2][1])', v)
                v = re.sub('z', '(self.baseVert[2] + xsi*self.dVerts[0][2] + eta*self.dVerts[1][2] + zet*self.dVerts[2][2])', v)
                self.funcXsiEtaZet = v
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




if __name__ == '__main__': 
    test0()
    test1()
    #test2()
    #test3()
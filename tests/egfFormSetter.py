import scipy
import re

class FormSetter:

    def __init__(self, order):
        """
        Constructor
        @param order order of the form (0, 1, 2, or 3)
        """
        self.order = order
        self.terms = {}
        self.formPat = re.compile(reduce(lambda x,y: x + '\^' + y, ['d\s*([xyz])\s*']*order))
        self.dVerts = []
        self.baseVert = None   

    def setExpression(self, expr):
        """
        Set form expression
        @param expr expression, eg x**2 * dx ^ dz + sin(pi*y)* dy ^ dx
        """
        coordName2Index = {'x': 0, 'y': 1, 'z': 2}
        # break the expresssion into independent terms
        continueFlag = True
        e = expr[:] # copy
        coords = []
        while continueFlag:
            m = re.search(self.formPat, e)
            if m:
                coords.append([coordName2Index[m.group(1 + o)] for o in self.order])
                e = re.sub(m.group(0), ',', 1) # replace the first occurrence only
            else:
                continueFlag = False
        self.terms = zip(coords, e.split(','))

    def setElement(self, *vertices):
        """
        Set the cell element
        @param vertices (ordert + 1)-tuple of vertices of the point, edge, face, or cell
        """
        assert len(vertices) == self.order + 1
        self.dVerts = [numpy.array(v) - numpy.array(vertices[0]) for v in vertices[1:]]
        self.baseVert = numpy.array(vertices[0])

    def evaluate(self):
        """
        Integrate the form over the element
        @return value
        """
        from numpy import sin, cos, tan, asin, acos, atan, atan2, pi, exp, log, log10, e
        res = 0
        if self.order == 0:
            for t in self.terms.values():
                res += eval(t)
        elif self.order == 1:
            for k, v in self.terms.items():
                k0, = k
                det = dVerts[0][k[0]]
                def f(xsi):
                    v = re.sub('x', '(self.baseVert[0] + xsi*self.dVerts[0][0])', v)
                    v = re.sub('y', '(self.baseVert[1] + xsi*self.dVerts[0][1])', v)
                    v = re.sub('z', '(self.baseVert[2] + xsi*self.dVerts[0][2])', v)
                    return eval(v)
                res += scipy.integrate.quad(f, 0., 1.)[0] * det
        elif self.order == 2:
            for k, v in self.terms.items():
                k0, k1 = k
                det = self.dVerts[0][k0] * self.dVerts[1][k1] - self.dVerts[0][k1] * self.dVerts[1][k0]
                def f(xsi, eta):
                    v = re.sub('x', '(self.baseVert[0] + xsi*self.dVerts[0][0] + eta*self.dVerts[1][0])', v)
                    v = re.sub('y', '(self.baseVert[1] + xsi*self.dVerts[0][1] + eta*self.dVerts[1][1])', v)
                    v = re.sub('z', '(self.baseVert[2] + xsi*self.dVerts[0][2] + eta*self.dVerts[1][2])', v)
                    return eval(v)
                def loEta(xsi):
                    return 0.0
                def hiEta(xsi):
                    return 1.0 - xsi
                res += scipy.integrate.dblquad(f, 0., 1., loEta, hiEta)[0] * det            
        elif self.order == 3:
            for k, v in self.terms.items():
                k0, k1, k2 = k
                det = scipy.linalg.det(jac)
                def f(xsi, eta, zet):
                    v = re.sub('x', '(self.baseVert[0] + xsi*self.dVerts[0][0] + eta*self.dVerts[1][0] + zet*self.dVerts[2][0])', v)
                    v = re.sub('y', '(self.baseVert[1] + xsi*self.dVerts[0][1] + eta*self.dVerts[1][1] + zet*self.dVerts[2][1])', v)
                    v = re.sub('z', '(self.baseVert[2] + xsi*self.dVerts[0][2] + eta*self.dVerts[1][2] + zet*self.dVerts[2][2])', v)
                    return eval(v)
                def loEta(xsi):
                    return 0.0
                def hiEta(xsi):
                    return 1.0 - xsi
                def loZet(xsi, eta):
                    return 0.0
                def hiZet(xsi, eta):
                    return 1.0 - xsi - eta             
                res += scipy.integrate.tplquad(f, 0., 1., loEta, hiEta, loZet, hiZet)[0] * det                        
        else:
            raise RunTimeError, 'ERROR: 0 <= order <= 3 but got {0}'.format(order)
        return res

###############################################################################
def test0():
    fs = FormSetter(0)
    fs.setExpression('x*y + z')
    fs.setElement([1., 2., 3.])
    res = fs.integrate()
    assert fabs(res - (1.*2. + 3.)) < 1.e-10

def test1():
    fs = FormSetter(0)
    fs.setExpression('x* dy + z**2 *dz')
    fs.setElement([1., 2., 3.], [1.1, 2.2, 3.3])
    res = fs.integrate()
    assert fabs(res - 0.21 - 2.979) < 1.e-10




if __name__ == '__main__': 
    test0()
    test1()
    #test2()
    #test3()
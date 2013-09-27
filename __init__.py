'''
Collider
======

The collider module contains classes which can be used to test membership
of a point in some space. See individual class documentation for details.
'''

__all__ = ('Collide2DPoly', )


class Collide2DPoly(object):
    ''' Collide2DPoly checks whether a point is within a polygon defined by a
    list of corner points.

    Based on http://alienryderflex.com/polygon/

    For example, a simple triangle::
    >>> collider = Collide2DPoly([10., 10., 20., 30., 30., 10.],\
                                 pre_compute=True)
    >>> (0.0, 0.0) in collider
    False
    >>> (20.0, 20.0) in collider
    True

    The constructor takes a list of x,y points in the form of [x1,y1,x2,y2...]
    as the points argument. These points define the corners of the
    polygon. The boundary is linearly interpolated between each set of points.
    The x, and y values must be floating points.
    The pre_compute argument, if True, will do some pre-computation to
    increase efficiency when comparing many points.
    '''

    def __init__(self, points, pre_compute=False, **kwargs):
        super(Collide2DPoly, self).__init__(**kwargs)
        self.points = points[:]
        if len(points) % 2:
            raise IndexError
        if len(points) < 4:
            self.points = None
            return
        self._min_x = min(points[0::2])
        self._max_x = max(points[0::2])
        self._min_y = min(points[1::2])
        self._max_y = max(points[1::2])

        self.pre_compute = pre_compute
        if not pre_compute:
            return
        count = len(points) / 2
        constant = [0, ] * count
        multiple = [0, ] * count
        j = count - 1
        for i in range(count):
            i_x = i * 2
            i_y = i_x + 1
            j_x = j * 2
            j_y = j_x + 1
            if points[j_y] == points[i_y]:
                constant[i] = points[i_x]
            else:
                constant[i] = (points[i_x] - points[i_y] * points[j_x] /
                               (points[j_y] - points[i_y]) +
                               points[i_y] * points[i_x] /
                               (points[j_y] - points[i_y]))
                multiple[i] = ((points[j_x] - points[i_x]) /
                               (points[j_y] - points[i_y]))
            j = i
        self._constant = constant
        self._multiple = multiple

    def collide_point(self, x, y):
        points = self.points
        if (not points) or not (self._min_x <= x <= self._max_x and
                                self._min_y <= y <= self._max_y):
            return False
        count = len(points) / 2
        j = count - 1
        odd = 0

        if self.pre_compute:
            constant = self._constant
            multiple = self._multiple
            for i in range(count):
                i_y = i * 2 + 1
                j_y = j * 2 + 1
                if (points[i_y] < y and points[j_y] >= y or
                    points[j_y] < y and points[i_y] >= y):
                    odd ^= y * multiple[i] + constant[i] < x
                j = i
        else:
            for i in range(count):
                i_x = i * 2
                i_y = i_x + 1
                j_x = j * 2
                j_y = j_x + 1
                if ((points[i_y] < y and points[j_y] >= y or
                     points[j_y] < y and points[i_y] >= y) and
                    (points[i_x] <= x or points[j_x] <= x)):
                    odd ^= (points[i_x] + (y - points[i_y])
                            / (points[j_y] - points[i_y])
                            * (points[j_x] - points[i_x]) < x)
                j = i
        return odd

    def __contains__(self, point):
        return self.collide_point(*point)


if __name__ == '__main__':
    from kivy.app import App
    from kivy.uix.widget import Widget
    from kivy.graphics import Line, Color, Rectangle
    from kivy.uix.button import Button
    from kivy.uix.boxlayout import BoxLayout
    from kivy.graphics.texture import Texture
    import itertools

    class CollideTester(Widget):

        def __init__(self, **kwargs):
            super(CollideTester, self).__init__(**kwargs)
            self.state = 'drawing'
            self.collider = None

        def on_touch_down(self, touch):
            if super(CollideTester, self).on_touch_down(touch):
                return True
            if not self.collide_point(*touch.pos):
                return False
            touch.grab(self)
            if self.state == 'drawing':
                with self.canvas:
                    Color(1, 0, 1, group='12345')
                    touch.ud['line'] = Line(points=[touch.x, touch.y],
                                            close=True, group='12345')

        def on_touch_move(self, touch):
            if touch.grab_current is not self:
                return super(CollideTester, self).on_touch_move(touch)
            if self.state == 'drawing':
                touch.ud['line'].points += [touch.x, touch.y]

        def on_touch_up(self, touch):
            if touch.grab_current is not self:
                return super(CollideTester, self).on_touch_up(touch)
            touch.ungrab(self)
            if self.state == 'drawing':
                self.state = 'testing'
                self.collider = Collide2DPoly(touch.ud['line'].points, True)
                collider = self.collider
                texture = Texture.create(size=self.size)
                inside = [255, 255, 0]
                outside = [255, 0, 255]
                x_off, y_off = self.pos
                width = int(self.width)
                height = int(self.height)
                buf = bytearray(width * height * 3)
                for x in range(width):
                    for y in range(height):
                        pos = (x + y * width) * 3
                        buf[pos:pos + 3] = (inside if (x + x_off, y + y_off)
                                            in collider else outside)
                texture.blit_buffer(bytes(buf), colorfmt='rgb',
                                    bufferfmt='ubyte')
                self.canvas.remove_group('12345')
                with self.canvas:
                    Rectangle(texture=texture, pos=self.pos, size=self.size,
                              group='12345')

    class TestApp(App):

        def build(self):
            box = BoxLayout(orientation='vertical')
            tester = CollideTester(size_hint_y=0.9)
            btn = Button(text='Clear', size_hint_y=0.1)
            box.add_widget(tester)
            box.add_widget(btn)

            def clear_state(*largs):
                tester.state = 'drawing'
                tester.collider = None
                tester.canvas.remove_group('12345')
            btn.bind(on_release=clear_state)
            return box

    TestApp().run()

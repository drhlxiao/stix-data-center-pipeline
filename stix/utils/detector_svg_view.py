#!/usr/bin/python3
#create a detector  svg file
from parser.utils import colour
import numpy as np

CANVAS_W = 1005
CANVAS_H = 1105
RFRAME = 450
RASPECT = 70
P0 = np.array([CANVAS_W / 2., CANVAS_H - RFRAME - 5])


def create_color_bar(x0, y0, width, height, max_value):
    path = (
        '<rect x="{x0}" y="{y0}" width="{width}" height="{height}"'
        'style="fill:rgb(250,250,250); stroke-width:0;stroke:rgb(0,0,255)" />'
    ).format(
        x0=x0, y=y0, width=width, height=height)
    num = len(colour.COLORS)
    dl = width / num
    for i, col in enumerate(colour.COLORS):
        x = dl * i + x0
        y = y0
        path += '<rect x="{}" y="{}" width="{}" height="{}"  style="fill:{}; stroke-width:1;stroke:{}" />'.format(
            x, y, dl, height, col, col)
    num_ticks = 10
    for i in range(0, num_ticks):
        x = dl * i + x0
        y = y0 - 20
        path += '<text x="{}" y="{}" > {} </text>'.format(
            x, y, max_value * i / num_ticks)

    return path


def create_STIX_svg(data=None):
    positions = [[135, 412.5], [135, 527.5], [135, 662.5], [135, 777.5],
                 [260, 297.5], [260, 412.5], [260, 527.5], [260, 662.5],
                 [260, 777.5], [260, 892.5], [385, 227.5], [385, 342.5],
                 [385, 457.5], [385, 732.5], [385, 847.5], [385, 962.5],
                 [510, 227.5], [510, 342.5], [510, 457.5], [510, 732.5],
                 [510, 847.5], [510, 962.5], [635, 297.5], [635, 412.5],
                 [635, 527.5], [635, 662.5], [635, 777.5], [635, 892.5],
                 [760, 412.5], [760, 527.5], [760, 662.5], [760, 777.5]]
    svg_start = '''
        <svg width="{cw}" height="{ch}">
      <circle
          style="opacity:0.1;fill:#0000ff;stroke:#0000ff;stroke-width:1.046;stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1"
         id="path5200"
         cx="{x0}"
         cy="{y0}"
         r="{r_frame}"/>
      <path  style="stroke-width:3;stroke:rgb(0,0,250)"  d="M{x0} {y_min} L{x0} {y_max} " />
      <path  style="stroke-width:3;stroke:rgb(0,0,250)"  d="M{x_min} {y0} L{x_max} {y0}" />
  <circle
     style="opacity:0.95999995;fill:#222b00;stroke:#0000ff;stroke-width:1.046;stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1"
     id="path5204"
     cx="{x0}"
     cy="{y0}"
     r="{r_aspect}" />
      '''.format(
        cw=CANVAS_W,
        ch=CANVAS_H,
        x0=P0[0],
        y0=P0[1],
        y_min=P0[1] - RFRAME,
        r_frame=RFRAME,
        r_aspect=RASPECT,
        y_max=P0[1] + RFRAME,
        x_min=P0[0] - RFRAME,
        x_max=P0[0] + RFRAME)

    svg_end = '   </svg>'

    color_bar = ''
    if data:
        color_bar = create_color_bar(10, 100, 1000, 30, max_value)
    template = svg_start + color_bar + '{}' + svg_end

    items = []
    for i, pos in enumerate(positions):
        pos[1] = pos[1] + 5

        items.append(create_one_detector_svg(i, pos, data, '{}'))

    return template.format('\n'.join(items))


def create_one_detector_svg(detector_id=0,
                            start=(0, 0),
                            data=[],
                            svg_template=''):
    """ data is a 32*12 array """
    max_value = 0
    if data:
        max_value = max(data)

    if not svg_template:
        svg_template = """
        <svg width="110" height="110">
        {}
        </svg>
        """
    group = ''

    pixel_template = '''
      <path id="{}" style="fill:{};stroke-width:1;stroke:rgb(200,200,200)"
      d="{}" />'''

    start_x = start[0]
    start_y = start[1]
    offset = np.array([start_x, start_y])
    fill_color = 'rgb(230,230,230)'

    pitch_x = 22
    pitch_y = 46

    big_p0_top = offset + np.array([6, 4])
    big_p0_bottom = offset + np.array([6, 4 + 92])

    big_pixel_top = 'h 22 v 46 h -11 v -4.5 h -11 Z'
    big_pixel_bottom = 'h 22 v -46 h -11 v 4.5 h -11 Z'

    small_p0 = offset + np.array([6, 50 - 4.5])
    small_pixel_path = 'h 11 v 9 h -11  Z'

    container = []
    #create big pixels

    guardring = '<rect x="{}" y="{}" width="100" height="100"  style="fill:rgb(255,255,255);stroke-width:1;stroke:rgb(0,0,0)" />'.format(
        offset[0], offset[1])
    x0 = offset[0]
    y0 = offset[1]
    #print("'det_{}':'M{} {} L {} {} L{} {} L{} {} Z ',".format(detector_id,x0,y0,x0+100,y0,x0+100,y0+100, x0,y0+100))
    container.append(guardring)
    text = '<text x="{}" y="{}" filled="red"> #{} </text>'.format(
        offset[0] + 40, offset[1] + 110, detector_id + 1)
    container.append(text)

    for i in range(0, 4):
        start = big_p0_top + np.array([i * pitch_x, 0])
        path = 'M {} {} {}'.format(start[0], start[1], big_pixel_top)
        pid = 'id-{}-{}'.format(detector_id, i)
        if data:
            fill_color = colour.get_color(data[detector_id * 12 + i],
                                          max_value)
        container.append(pixel_template.format(pid, fill_color, path))

    for i in range(0, 4):
        start = big_p0_bottom + np.array([i * pitch_x, 0])
        path = 'M {} {} {}'.format(start[0], start[1], big_pixel_bottom)
        pid = 'id-{}-{}'.format(detector_id, i + 4)
        if data:
            fill_color = colour.get_color(data[detector_id * 12 + i],
                                          max_value)
        container.append(pixel_template.format(pid, fill_color, path))

    for i in range(0, 4):
        start = small_p0 + np.array([i * pitch_x, 0])
        path = 'M {} {} {}'.format(start[0], start[1], small_pixel_path)
        pid = 'id-{}-{}'.format(detector_id, i + 8)
        if data:
            fill_color = colour.get_color(data[detector_id * 12 + i],
                                          max_value)
        container.append(pixel_template.format(pid, fill_color, path))
    group = '<g> {} </g>'.format('\n'.join(container))
    return svg_template.format(group)


if __name__ == '__main__':
    with open('detector.svg', 'w') as f:
        f.write(create_STIX_svg())

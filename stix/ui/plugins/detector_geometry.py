import colour
import numpy as np


def create_STIX_svg(data=None, filename=''):
    positions = [(40.000, 270.000), (40.000, 390.000), (40.000, 550.000),
                 (40.000, 670.000), (180.000, 150.000), (180.000, 270.000),
                 (180.000, 390.000), (180.000, 550.000), (180.000, 670.000),
                 (180.000, 790.000), (310.000, 70.000), (310.000, 190.000),
                 (310.000, 310.000), (310.000, 630.000), (310.000, 750.000),
                 (310.000, 870.000), (430.000, 70.000), (430.000, 190.000),
                 (430.000, 310.000), (430.000, 630.000), (430.000, 750.000),
                 (430.000, 870.000), (570.000, 150.000), (570.000, 270.000),
                 (570.000, 390.000), (570.000, 550.000), (570.000, 670.000),
                 (570.000, 790.000), (710.000, 270.000), (710.000, 390.000),
                 (710.000, 550.000), (710.000, 670.000)]
    template = '''
        <svg width="861" height="1020">
    <rect x="10" y="10" width="840" height="1010"  style="fill:rgb(255,255,255); stroke-width:3;stroke:rgb(0,0,255)" />
      <path  style="stroke-width:3;stroke:rgb(0,0,250)"  d="M420 10 L420 1020" />
      <path  style="stroke-width:3;stroke:rgb(0,0,250)"  d="M10 525 L851 525" />
        {}
        </svg>'''
    items = []
    for i, pos in enumerate(positions):
        detector_data = None
        if data:
            detector_data = data[i]

        items.append(create_one_detector_svg(i, pos, detector_data, '{}'))
    svg=template.format('\n'.join(items))

    if filename:
        with open(filename,'w') as f:
            f.write(svg)

    return svg



def create_one_detector_svg(detector_id=0,
                            start=(0, 0),
                            data=None,
                            svg_template=''):
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
    container.append(guardring)
    text = '<text x="{}" y="{}" filled="red"> #{} </text>'.format(
        offset[0] + 40, offset[1] + 110, detector_id + 1)
    container.append(text)

    colors = colour.get_colors(data)

    for i in range(0, 4):
        start = big_p0_top + np.array([i * pitch_x, 0])
        path = 'M {} {} {}'.format(start[0], start[1], big_pixel_top)
        pid = 'id-{}-{}'.format(detector_id, i)
        if colors:
            fill_color = colors[i]
        container.append(pixel_template.format(pid, fill_color, path))

    for i in range(0, 4):
        start = big_p0_bottom + np.array([i * pitch_x, 0])
        path = 'M {} {} {}'.format(start[0], start[1], big_pixel_bottom)
        pid = 'id-{}-{}'.format(detector_id, i + 4)
        if colors:
            fill_color = colors[i + 4]
        container.append(pixel_template.format(pid, fill_color, path))

    for i in range(0, 4):
        start = small_p0 + np.array([i * pitch_x, 0])
        path = 'M {} {} {}'.format(start[0], start[1], small_pixel_path)
        pid = 'id-{}-{}'.format(detector_id, i + 8)
        if colors:
            fill_color = colors[i + 8]

        container.append(pixel_template.format(pid, fill_color, path))
    group = '<g> {} </g>'.format('\n'.join(container))
    return svg_template.format(group)


if __name__ == '__main__':
    with open('detector.svg', 'w') as f:
        f.write(create_STIX_svg())
        ##create_one_detector_svg(
        #    data=[23, 5, 43, 543, 53, 33, 43, 33, 34, 54, 55, 44]))

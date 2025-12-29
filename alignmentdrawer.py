import turtle
import tkinter as tk
from PIL import Image

CONFIG = {
    'fastafile': 'seq.fasta',
    'referenceseqid': 'refseq',
    'outfile': 'outfile.tiff',
    'line_spacing': 30,          
    'base_line_width': 15,       
    'gap_line_width': 5,         
    'left_margin': 180,          
    'right_margin': 50,          
    'font_spec': ('Arial', 16),  
    'dpi': 600,                  
    
    'colors': {
        'match': (100, 200, 240),        
        'mismatch': (255, 180, 200),   
        'gap': (0, 0, 0),         
        'background': (255, 255, 255) 
    },
    
    'ruler': {
        'enabled': True,         
        'height': 10,            
        'tick_interval': 100,    
        'minor_tick': 50,        
        'line_width': 2,         
        'font': ('Arial', 16),   
        'color': (0, 0, 0),        
        'tick_length': 5        
    }
}

def get_color(color_def):
    if isinstance(color_def, str) and color_def.startswith('colors.'):
        color_name = color_def.split('.')[1]
        rgb = CONFIG['colors'][color_name]
        return "#%02x%02x%02x" % rgb  
    elif isinstance(color_def, (tuple, list)):
        return "#%02x%02x%02x" % tuple(color_def)
    return color_def  

def calculate_canvas_size(fasta_dict, seq_length):
    num_sequences = len(fasta_dict)
    height = num_sequences * CONFIG['line_spacing'] + 100
    
    
    if CONFIG['ruler']['enabled']:
        height += CONFIG['line_spacing']  
    
    width = seq_length + CONFIG['left_margin'] + CONFIG['right_margin']
    return width, height

def read_fasta(fastafile):
    fasta_dict = {}
    with open(fastafile) as f:
        for record in f.read().split('>')[1:]:
            header, *seq_lines = record.split('\n')
            seqid = header.split()[0]
            sequence = ''.join(seq_lines).upper()
            fasta_dict[seqid] = sequence
    return fasta_dict

def draw_full_ruler(canvas, canvas_width, canvas_height, seq_length, last_seq_y):
    if not CONFIG['ruler']['enabled']:
        return
    
    ruler_config = CONFIG['ruler']
    start_x = -canvas_width//2 + CONFIG['left_margin']
    end_x = start_x + seq_length
    
    
    ruler_y = last_seq_y + CONFIG['line_spacing']
    
    
    drawn_positions = set()
    
    
    for pos in range(0, seq_length + 1, ruler_config['tick_interval']):
        tick_x = start_x + pos
        drawn_positions.add(pos)
        
        canvas.create_line(
            tick_x, ruler_y - ruler_config['tick_length'],
            tick_x, ruler_y + ruler_config['tick_length'],
            width=ruler_config['line_width'],
            fill=get_color(ruler_config['color'])
        )
        
        canvas.create_text(
            tick_x, ruler_y + ruler_config['tick_length'] + 5,
            text=str(pos),
            anchor=tk.N,
            font=ruler_config['font'],
            fill=get_color(ruler_config['color'])
        )
    
    
    if seq_length not in drawn_positions:
        tick_x = end_x
        
        canvas.create_line(
            tick_x, ruler_y - ruler_config['tick_length'],
            tick_x, ruler_y + ruler_config['tick_length'],
            width=ruler_config['line_width'],
            fill=get_color(ruler_config['color'])
        )
        
        canvas.create_text(
            tick_x ,  
            ruler_y,
            text=str(seq_length),
            anchor=tk.W,
            font=ruler_config['font'],
            fill=get_color(ruler_config['color'])
        )
    
    
    for pos in range(0, seq_length + 1, ruler_config['minor_tick']):
        if pos % ruler_config['tick_interval'] != 0 and pos != seq_length:
            tick_x = start_x + pos
            canvas.create_line(
                tick_x, ruler_y - ruler_config['tick_length']//2,
                tick_x, ruler_y + ruler_config['tick_length']//2,
                width=ruler_config['line_width'],
                fill=get_color(ruler_config['color'])
            )
    
    
    canvas.create_line(
        start_x, ruler_y,
        end_x, ruler_y,
        width=ruler_config['line_width'],
        fill=get_color(ruler_config['color'])
    )


def main():
    
    fastadict = read_fasta(CONFIG['fastafile'])
    refseq = fastadict[CONFIG['referenceseqid']]
    seq_length = len(refseq)
    
    
    canvas_width, canvas_height = calculate_canvas_size(fastadict, seq_length)
    
    
    screen = turtle.Screen()
    screen.setup(width=canvas_width, height=canvas_height)
    turtle.tracer(0)
    turtle.ht()
    canvas = screen.getcanvas()
    
    
    start_x = -canvas_width//2 + CONFIG['left_margin']
    start_y = -canvas_height//2 + 50  
    
    
    last_seq_y = None  
    for seqid, sequence in fastadict.items():
        current_x = start_x
        current_y = start_y
        
        
        canvas.create_text(
            -canvas_width//2 + 10, current_y,
            text=seqid,
            anchor=tk.W,
            font=CONFIG['font_spec'],
            fill='black'
        )
        
        
        bi = 0
        for base in sequence:
            if base == '-':
                color = get_color(CONFIG['colors']['gap'])
                width = CONFIG['gap_line_width']
            elif base == refseq[bi]:
                color = get_color(CONFIG['colors']['match'])
                width = CONFIG['base_line_width']
            else:
                color = get_color(CONFIG['colors']['mismatch'])
                width = CONFIG['base_line_width']
            
            canvas.create_line(
                current_x, current_y,
                current_x + 1, current_y,
                width=width,
                fill=color,
                capstyle=tk.BUTT,
                joinstyle=tk.MITER
            )
            current_x += 1
            bi += 1
        
        start_y += CONFIG['line_spacing']  
        last_seq_y = current_y  

    
    
    draw_full_ruler(canvas, canvas_width, canvas_height, seq_length, last_seq_y)
    
    
    def export_tiff():
        ps_file = "temp.ps"
        canvas.postscript(
            file=ps_file,
            width=canvas_width,
            height=canvas_height,
            x=-canvas_width//2,
            y=-canvas_height//2
        )
        img = Image.open(ps_file)
        img.save(CONFIG['outfile'], "TIFF", dpi=(CONFIG['dpi'], CONFIG['dpi']))

    
    turtle.update()
    export_tiff()


if __name__ == '__main__':
    main()
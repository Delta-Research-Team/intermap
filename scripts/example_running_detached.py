# Created by rglez at 7/18/25
from bitarray.util import sc_decode

from intermap import runner

cfg_path = '/media/rglez/Roy2TB/Dropbox/RoyData/intermap/tutorial-mayank/channel-5.cfg'
bit_dict = runner.execute(cfg_path=cfg_path)
bit_dict = {k: sc_decode(v) for k, v in bit_dict.items()}

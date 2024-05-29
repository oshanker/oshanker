import base_transformer_shanker.data
from torch.utils.data import DataLoader
from torch.nn.utils.rnn import pad_sequence
from base_transformer_shanker.data import GenerateDataset as GD
from base_transformer_shanker.functions import gd_collate_fn, tokens_to_str
from base_transformer_shanker.constants import PAD_IDX, SOS_IDX, EOS_IDX


def runpersonMain():
	train_iter = GD(6, pad_idx=PAD_IDX, sos_idx=SOS_IDX, eos_idx=EOS_IDX)
	dataloader_train = DataLoader(train_iter, batch_size=3, collate_fn=gd_collate_fn)
	s, t = next(iter(dataloader_train))
	print(s[:, ...])
	print(t[:, ...])
	print("s.size", s.size())
	print("?", tokens_to_str(s[0, :]))
	print("?", tokens_to_str(t[0, :]))
	
	for s, t in iter(dataloader_train):
	     print(s.size())


print(dir(base_transformer_shanker.data))

runpersonMain()


# export PYTHONPATH="$PWD"
# python3 run.py
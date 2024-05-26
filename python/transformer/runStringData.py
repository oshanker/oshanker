import base_transformer_shanker.data
from torch.utils.data import DataLoader


PAD_IDX = base_transformer_shanker.data.ReverseDataset.PAD_IDX
SOS_IDX = base_transformer_shanker.data.ReverseDataset.SOS_IDX
EOS_IDX = base_transformer_shanker.data.ReverseDataset.EOS_IDX

def tokens_to_str(tokens):
	return "".join([chr(x+94) for x in tokens])


def runpersonMain():
	train_iter = base_transformer_shanker.data.ReverseDataset(6, pad_idx=PAD_IDX, sos_idx=SOS_IDX, eos_idx=EOS_IDX)
	dataloader_train = DataLoader(train_iter, batch_size=3)
	s, t = next(iter(dataloader_train))
	print(s[:, ...])
	print(t[:, ...])
	print("s.size", s.size())
	print("?", tokens_to_str(s[0, :]))
	print("?", tokens_to_str(t[0, :]))


print(dir(base_transformer_shanker.data))

runpersonMain()


# export PYTHONPATH="$PWD"
# python3 run.py
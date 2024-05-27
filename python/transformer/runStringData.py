import base_transformer_shanker.data
from torch.utils.data import DataLoader
from torch.nn.utils.rnn import pad_sequence
from base_transformer_shanker.data import GenerateDataset as GD

PAD_IDX = GD.PAD_IDX
SOS_IDX = GD.SOS_IDX
EOS_IDX = GD.EOS_IDX

def tokens_to_str(tokens):
	return "".join([chr(x+94) for x in tokens])

def collate_fn(batch):
    """ 
    This function pads inputs with PAD_IDX to have batches of equal length
    """
    src_batch, tgt_batch = [], []
    for src_sample, tgt_sample in batch:
        src_batch.append(src_sample)
        tgt_batch.append(tgt_sample)

    src_batch = pad_sequence(src_batch, padding_value=PAD_IDX, batch_first=True)
    tgt_batch = pad_sequence(tgt_batch, padding_value=PAD_IDX, batch_first=True)
    return src_batch, tgt_batch


def runpersonMain():
	train_iter = GD(6, pad_idx=PAD_IDX, sos_idx=SOS_IDX, eos_idx=EOS_IDX)
	dataloader_train = DataLoader(train_iter, batch_size=3, collate_fn=collate_fn)
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
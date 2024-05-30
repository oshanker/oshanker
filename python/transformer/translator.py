import torch
import torch.nn as nn
from transformer import Transformer
from base_transformer_shanker.constants import args, SOS_IDX, EOS_IDX

#from base_transformer_shanker.data import GenerateDataset as GD


class Translator(nn.Module):
    """
    https://pytorch.org/tutorials/beginner/saving_loading_models.html
    """

    def __init__(self, transformer):
        super(Translator, self).__init__()
        self.transformer = transformer
    
    @staticmethod
    def str_to_tokens(s):
        return [ord(z)-97+3 for z in s]
    
    @staticmethod
    def tokens_to_str(tokens):
        return "".join([chr(x+94) for x in tokens])
    
    def __call__(self, sentence, max_length=None, pad=False):
        
        x = torch.tensor(self.str_to_tokens(sentence))
        x = torch.cat([torch.tensor([SOS_IDX]), x, torch.tensor([EOS_IDX])]).unsqueeze(0)
        
        encoder_output, mask = self.transformer.encode(x) # (B, S, E)
        
        if not max_length:
            max_length = x.size(1)
            
        outputs = torch.ones((x.size()[0], max_length)).type_as(x).long() * SOS_IDX
        
        for step in range(1, max_length):
            y = outputs[:, :step]
            probs = self.transformer.decode(y, encoder_output)
            output = torch.argmax(probs, dim=-1)
            currentString = Translator.tokens_to_str(y[0])
            print(f"Knowing {currentString} we output {output[:, -1]}")
            if output[:, -1].detach().numpy() in (EOS_IDX, SOS_IDX):
                break
            outputs[:, step] = output[:, -1]
            
        
        return self.tokens_to_str(outputs[0])

def runperson():
    # Define model here
    model = Transformer(**args)
    path = "../out/forward.pt"
    model.load_state_dict(torch.load(path))
    model.eval()
    translator = Translator(model)
    sentence = "helloworld"
    print(sentence)
    out = translator(sentence)
    print(out)
    
if __name__ == "__main__":
    runperson()

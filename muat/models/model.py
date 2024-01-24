import math

from loguru import logger
import torch
import torch.nn as nn
from torch.nn import functional as F
import pdb

class ModelConfig:
    """ base GPT config, params common to all GPT versions """
    embd_pdrop = 0.1
    resid_pdrop = 0.1
    attn_pdrop = 0.1

    def __init__(self,vocab_size, block_size, num_class,position_size = 1, ges_size=1,dnn_input=1,epi_size=1,args=None, **kwargs):
        self.vocab_size = vocab_size
        self.block_size = block_size
        self.num_class = num_class
        self.position_size = position_size
        self.ges_size = ges_size
        self.args = args
        self.dnn_input=dnn_input
        self.epi_size=epi_size

        for k,v in kwargs.items():
            setattr(self, k, v)

class Block(nn.Module):
    """ an unassuming Transformer block """

    def __init__(self, config):
        super().__init__()
        self.ln1 = nn.LayerNorm(config.n_embd)
        self.ln2 = nn.LayerNorm(config.n_embd)
        self.attn = CausalSelfAttention(config) #TODO: https://pytorch.org/tutorials/intermediate/scaled_dot_product_attention_tutorial.html#causal-self-attention
        self.mlp = nn.Sequential(
            nn.Linear(config.n_embd, 4 * config.n_embd),
            nn.GELU(),
            nn.Linear(4 * config.n_embd, config.n_embd),
            nn.Dropout(config.resid_pdrop),
        )

    def forward(self, x):
        x = x + self.attn(self.ln1(x))
        x = x + self.mlp(self.ln2(x))
        return x

class DNN_GX(nn.Module):
    """
    Transformer for classifying sequences
    """
    def __init__(self, config):
        super().__init__()

        self.dnn = nn.Sequential(nn.Linear(config.dnn_input, 10000),
                                   nn.ReLU(),
                                   nn.Linear(10000, 2500),
                                   nn.ReLU(),
                                   nn.Linear(2500, 500),
                                   nn.ReLU(),
                                   nn.Linear(500,24),
                                   nn.ReLU(),
                                   )

        self.toprobs = nn.Linear(24, config.num_class)

    def forward(self, x,targets=None,vis=None):

        x = x[0]
        x = x[0]


        features = self.dnn(x)

        if vis:
            return features

        logits = self.toprobs(features)

        # if we are given some desired targets also calculate the loss
        loss = None
        if targets is not None:
            loss = F.cross_entropy(logits.view(-1, logits.size(-1)), targets.view(-1))

        return logits, loss


class EmbFC(nn.Module):
    """
    Transformer for classifying sequences
    """
    def __init__(self, config):
        super().__init__()

        self.num_tokens, self.max_pool = config.vocab_size, False

        self.token_embedding = nn.Embedding(config.vocab_size, config.n_embd)

        self.fcemb = nn.Sequential(nn.Linear(config.block_size, config.n_embd),
                                   nn.ReLU(),
                                   nn.Linear(config.n_embd, config.n_embd),
                                   nn.ReLU(),
                                   nn.Linear(config.n_embd, 1),
                                   nn.ReLU(),
                                   )


        self.toprobs = nn.Linear(config.n_embd, config.num_class)

        self.do = nn.Dropout(config.embd_pdrop)

    def forward(self, x,targets=None,vis=None):

        x = x[0]

        tokens = self.token_embedding(x)

        tokens = tokens.permute(0, 2, 1)
        tokens = tokens.squeeze()

        tokens = self.fcemb(tokens)
 


        tokens = tokens.permute(1, 0)

        logits = self.toprobs(tokens)

        # if we are given some desired targets also calculate the loss
        loss = None
        if targets is not None:
            loss = F.cross_entropy(logits.view(-1, logits.size(-1)), targets.view(-1))

        return logits, loss

class EmbFCMean(nn.Module):
    """
    Transformer for classifying sequences
    """
    def __init__(self, config):
        super().__init__()
        self.num_tokens, self.max_pool = config.vocab_size, False
        self.token_embedding = nn.Embedding(config.vocab_size, config.n_embd)
        self.fcemb = nn.Sequential(
            nn.Linear(config.block_size, config.n_embd),
            nn.ReLU(),
            nn.Linear(config.n_embd, config.n_embd),
            nn.ReLU(),
            nn.Linear(config.n_embd, 1),
            nn.ReLU()
        )
        self.toprobs = nn.Linear(config.n_embd, config.num_class)
        self.do = nn.Dropout(config.embd_pdrop)

    def forward(self, x,targets=None,vis=None):

        x = x[0]
        #pdb.set_trace()

        tokens = self.token_embedding(x)

        tokens = tokens.permute(0, 2, 1)
        tokens = tokens.squeeze()

        tokens = self.fcemb(tokens)
 


        tokens = tokens.permute(1, 0)

        logits = self.toprobs(tokens)

        # if we are given some desired targets also calculate the loss
        loss = None
        if targets is not None:
            loss = F.cross_entropy(logits.view(-1, logits.size(-1)), targets.view(-1))

        return logits, loss

class EmbFCPos(nn.Module):
    """
    Transformer for classifying sequences
    """
    def __init__(self, config):
        super().__init__()

        self.num_tokens, self.max_pool = config.vocab_size, False

        self.token_embedding = nn.Embedding(config.vocab_size, config.n_embd)
        self.position_embedding = nn.Embedding(config.position_size, config.n_embd,padding_idx=0)

        self.fcemb = nn.Sequential(nn.Linear(config.block_size, config.n_embd + config.n_embd),
                                   nn.ReLU(),
                                   nn.Linear(config.n_embd + config.n_embd, config.n_embd),
                                   nn.ReLU(),
                                   nn.Linear(config.n_embd, 1),
                                   nn.ReLU(),
                                   )


        self.toprobs = nn.Linear(config.n_embd, config.num_class)

        self.do = nn.Dropout(config.embd_pdrop)

    def forward(self, x,targets=None,vis=None):

        triplet = x[0]
        #pdb.set_trace()
        pos = x[1]

        #pdb.set_trace()

        tokens = self.token_embedding(triplet)
        position = self.position_embedding(pos)

        x = torch.cat((tokens,position),axis=2)

        tokens = tokens.permute(0, 2, 1)
        tokens = tokens.squeeze()

        tokens = self.fcemb(tokens)

        tokens = tokens.permute(1, 0)

        logits = self.toprobs(tokens)

        # if we are given some desired targets also calculate the loss
        loss = None
        if targets is not None:
            loss = F.cross_entropy(logits.view(-1, logits.size(-1)), targets.view(-1))

        return logits, loss

class EmbFCPosGES(nn.Module):
    """
    Transformer for classifying sequences
    """
    def __init__(self, config):
        super().__init__()

        self.num_tokens, self.max_pool = config.vocab_size, False

        self.token_embedding = nn.Embedding(config.vocab_size, config.n_embd)
        self.position_embedding = nn.Embedding(config.position_size, config.n_embd,padding_idx=0)
        self.ges_embedding = nn.Embedding(config.ges_size+1, 4,padding_idx=0)

        self.fcemb = nn.Sequential(nn.Linear(config.block_size, config.n_embd + config.n_embd + 4),
                                   nn.ReLU(),
                                   nn.Linear(config.n_embd + config.n_embd + 4, config.n_embd),
                                   nn.ReLU(),
                                   nn.Linear(config.n_embd, 1),
                                   nn.ReLU(),
                                   )


        self.toprobs = nn.Linear(config.n_embd, config.num_class)

        self.do = nn.Dropout(config.embd_pdrop)

    def forward(self, x,targets=None,vis=None):

        triplet = x[0]
        pos = x[1]
        ges = x[2]
        #pdb.set_trace()

        tokens = self.token_embedding(triplet)
        position = self.position_embedding(pos)
        gesemb = self.ges_embedding(ges)

        x = torch.cat((tokens,position,gesemb),axis=2)

        tokens = tokens.permute(0, 2, 1)
        tokens = tokens.squeeze()

        tokens = self.fcemb(tokens)

        tokens = tokens.permute(1, 0)

        logits = self.toprobs(tokens)

        # if we are given some desired targets also calculate the loss
        loss = None
        if targets is not None:
            loss = F.cross_entropy(logits.view(-1, logits.size(-1)), targets.view(-1))

        return logits, loss


class MuAtMotif(nn.Module):
    """
    Transformer for classifying sequences
    """
    def __init__(self, config):
        super().__init__()

        self.num_tokens, self.max_pool = config.vocab_size, False

        self.token_embedding = nn.Embedding(config.vocab_size, config.n_embd,padding_idx=0)

        tblocks = []
        for i in range(config.n_layer):
            tblocks.append(
                TransformerBlock(emb=config.n_embd, heads=config.n_head, seq_length=config.block_size, mask=False, dropout=config.attn_pdrop))

        self.tblocks = nn.Sequential(*tblocks)

        self.toprobs = nn.Linear(config.n_embd, config.num_class)

        self.do = nn.Dropout(config.embd_pdrop)

    def forward(self, x,targets=None,vis=None):

        x = x[0]
        #pdb.set_trace()
    
        tokens = self.token_embedding(x)

        b, t, e = tokens.size()
        
        x = self.do(tokens)

        x = self.tblocks(x)

        if vis:
            return x

        x = x.max(dim=1)[0] if self.max_pool else x.mean(dim=1) # pool over the time dimension

        logits = self.toprobs(x)

        # if we are given some desired targets also calculate the loss
        loss = None
        if targets is not None:
            loss = F.cross_entropy(logits.view(-1, logits.size(-1)), targets.view(-1))
        return logits, loss

class MuAtEpiF(nn.Module):
    """
    Transformer for classifying sequences
    """
    def __init__(self, config):
        super().__init__()

        self.num_tokens, self.max_pool = config.vocab_size, False

        self.device = 'cpu'
        if torch.cuda.is_available():
            self.device = torch.cuda.current_device()

        self.epiblocks = nn.Embedding(config.epi_size, config.epi_emb,padding_idx=0)
        self.epi_emb = config.epi_emb

        tblocks = []
        for i in range(config.n_layer):
            tblocks.append(
                TransformerBlock(emb= int(config.epi_emb * 145), heads=config.n_head, seq_length=config.block_size, mask=False, dropout=config.attn_pdrop))

        self.tblocks = nn.Sequential(*tblocks)

        self.tofeature = nn.Sequential(nn.Linear(int(config.epi_emb * 145), 24),
                                        nn.ReLU())

        self.toprobs = nn.Linear(24, config.num_class)

        self.do = nn.Dropout(config.embd_pdrop)

    def forward(self, x,targets=None,vis=None):

        x = x[0]

        batch = x.shape[0]

        all_emb = torch.zeros(batch,x.shape[1], self.epi_emb*x.shape[2]).to(self.device)
        start = 0
        end = self.epi_emb
        for i in range(x.shape[2]):
            oneinput = x[:,:,i]
            oneres = self.epiblocks.to(self.device)(oneinput)
            all_emb[:,:,start:end] = oneres
            start = start + self.epi_emb
            end = end + self.epi_emb

        x = self.tblocks(all_emb)

        x = x.max(dim=1)[0] if self.max_pool else x.mean(dim=1) # pool over the time dimension

        feature = self.tofeature(x)

        if vis:
            return feature
        else:
            logits = self.toprobs(feature)

        # if we are given some desired targets also calculate the loss
        loss = None
        if targets is not None:
            loss = F.cross_entropy(logits.view(-1, logits.size(-1)), targets.view(-1))

        return logits, loss

class MuAtMotifPositionGESEpiF(nn.Module):
    """
    Transformer for classifying sequences
    """
    def __init__(self, config):
        super().__init__()

        self.device = 'cpu'
        if torch.cuda.is_available():
            self.device = torch.cuda.current_device()

        self.num_tokens, self.max_pool = config.vocab_size, False

        self.token_embedding = nn.Embedding(config.vocab_size, config.n_embd,padding_idx=0)
        self.position_embedding = nn.Embedding(config.position_size, config.n_embd,padding_idx=0)
        self.ges_embedding = nn.Embedding(config.ges_size+1, 4,padding_idx=0)


        self.epiblocks = []
        for i in range(145):
            self.epiblocks.append(nn.Embedding(19, 2,padding_idx=0))

        tblocks = []
        for i in range(config.n_layer):
            tblocks.append(
                TransformerBlock(emb=int(config.n_embd + config.n_embd + 4 + (2 * 145)), heads=config.n_head, seq_length=config.block_size, mask=False, dropout=config.attn_pdrop))

        self.tblocks = nn.Sequential(*tblocks)

        self.tofeature = nn.Sequential(nn.Linear(int(config.n_embd + config.n_embd + 4 + (2 * 145)), 24),
                                        nn.ReLU())

        self.toprobs = nn.Linear(24, config.num_class)

        self.do = nn.Dropout(config.embd_pdrop)

    def forward(self, x,targets=None,vis=None):

        triplettoken = x[0]
        postoken = x[1]
        gestoken = x[2]

        epitoken = x[3]

        batch = triplettoken.shape[0]

        all_emb = torch.zeros(batch,epitoken.shape[1], 2*epitoken.shape[2]).to(self.device)
        
        start = 0
        end = 2
        for i in range(epitoken.shape[2]):
            oneinput = epitoken[:,:,i]
            oneres = self.epiblocks[i].to(self.device)(oneinput)
            all_emb[:,:,start:end] = oneres
            start = start + 2
            end = end + 2
    
        tokens = self.token_embedding(triplettoken)
        positions = self.position_embedding(postoken)
        ges = self.ges_embedding(gestoken)

        x = torch.cat((tokens,positions,ges,all_emb),axis=2)

        x = self.do(x)
        #pdb.set_trace()
        
        x = self.tblocks(x)

        x = x.max(dim=1)[0] if self.max_pool else x.mean(dim=1) # pool over the time dimension

        feature = self.tofeature(x)

        if vis:
            return feature
        else:
            logits = self.toprobs(feature)

        # if we are given some desired targets also calculate the loss
        loss = None
        if targets is not None:
            loss = F.cross_entropy(logits.view(-1, logits.size(-1)), targets.view(-1))

        return logits, loss

class MuAtMotifPositionGESEpiF_OneEpiEmb(nn.Module):
    """
    Transformer for classifying sequences
    """
    def __init__(self, config):
        super().__init__()

        self.device = 'cpu'
        if torch.cuda.is_available():
            self.device = torch.cuda.current_device()

        self.num_tokens, self.max_pool = config.vocab_size, False

        self.token_embedding = nn.Embedding(config.vocab_size, config.n_embd,padding_idx=0)
        self.position_embedding = nn.Embedding(config.position_size, config.n_embd,padding_idx=0)
        self.ges_embedding = nn.Embedding(config.ges_size+1, 4,padding_idx=0)
        #pdb.set_trace()
        self.epiblocks = nn.Embedding(config.epi_size, config.epi_emb,padding_idx=0)

        self.epi_emb = config.epi_emb
        
        tblocks = []
        for i in range(config.n_layer):
            tblocks.append(
                TransformerBlock(emb=int(config.n_embd + config.n_embd + 4 + (config.epi_emb * 145)), heads=config.n_head, seq_length=config.block_size, mask=False, dropout=config.attn_pdrop))

        self.tblocks = nn.Sequential(*tblocks)

        self.tofeature = nn.Sequential(nn.Linear(int(config.n_embd + config.n_embd + 4 + (config.epi_emb * 145)), 24),
                                        nn.ReLU())

        self.toprobs = nn.Linear(24, config.num_class)

        self.do = nn.Dropout(config.embd_pdrop)

    def forward(self, x,targets=None,vis=None):

        triplettoken = x[0]
        postoken = x[1]
        gestoken = x[2]

        epitoken = x[3]

        batch = triplettoken.shape[0]

        all_emb = torch.zeros(batch,epitoken.shape[1], self.epi_emb*epitoken.shape[2]).to(self.device)
        
        start = 0
        end = self.epi_emb
        for i in range(epitoken.shape[2]):
            oneinput = epitoken[:,:,i]
            oneres = self.epiblocks.to(self.device)(oneinput)
            all_emb[:,:,start:end] = oneres
            start = start + self.epi_emb
            end = end + self.epi_emb
    
        tokens = self.token_embedding(triplettoken)
        positions = self.position_embedding(postoken)
        ges = self.ges_embedding(gestoken)

        x = torch.cat((tokens,positions,ges,all_emb),axis=2)

        x = self.do(x)
        #pdb.set_trace()
        
        x = self.tblocks(x)

        x = x.max(dim=1)[0] if self.max_pool else x.mean(dim=1) # pool over the time dimension

        feature = self.tofeature(x)

        if vis:
            return feature
        else:
            logits = self.toprobs(feature)

        # if we are given some desired targets also calculate the loss
        loss = None
        if targets is not None:
            loss = F.cross_entropy(logits.view(-1, logits.size(-1)), targets.view(-1))

        return logits, loss


class MuAtMotifPositionGESEpiF_EpiCompress(nn.Module):
    """
    Transformer for classifying sequences
    """
    def __init__(self, config):
        super().__init__()

        self.device = 'cpu'
        if torch.cuda.is_available():
            self.device = torch.cuda.current_device()

        self.num_tokens, self.max_pool = config.vocab_size, False

        self.token_embedding = nn.Embedding(config.vocab_size, config.n_embd,padding_idx=0)
        self.position_embedding = nn.Embedding(config.position_size, config.n_embd,padding_idx=0)
        self.ges_embedding = nn.Embedding(config.ges_size+1, 4,padding_idx=0)
        #pdb.set_trace()
        self.epiblocks = nn.Embedding(config.epi_size, config.epi_emb,padding_idx=0)
        self.epi_emb=config.epi_emb

        self.epi_enc = nn.Linear((config.epi_emb * 145), config.n_embd)
        
        tblocks = []
        for i in range(config.n_layer):
            tblocks.append(
                TransformerBlock(emb=int(config.n_embd + config.n_embd + 4 + (config.n_embd)), heads=config.n_head, seq_length=config.block_size, mask=False, dropout=config.attn_pdrop))

        self.tblocks = nn.Sequential(*tblocks)

        self.tofeature = nn.Sequential(nn.Linear(int(config.n_embd + config.n_embd + 4 + (config.n_embd)), 24),
                                        nn.ReLU())

        self.toprobs = nn.Linear(24, config.num_class)

        self.do = nn.Dropout(config.embd_pdrop)

    def forward(self, x,targets=None,vis=None):

        triplettoken = x[0]
        postoken = x[1]
        gestoken = x[2]

        epitoken = x[3]

        batch = triplettoken.shape[0]

        all_emb = torch.zeros(batch,epitoken.shape[1], self.epi_emb*epitoken.shape[2]).to(self.device)
        
        start = 0
        end = self.epi_emb
        for i in range(epitoken.shape[2]):
            oneinput = epitoken[:,:,i]
            oneres = self.epiblocks.to(self.device)(oneinput)
            all_emb[:,:,start:end] = oneres
            start = start + self.epi_emb
            end = end + self.epi_emb
    
        tokens = self.token_embedding(triplettoken)
        positions = self.position_embedding(postoken)
        ges = self.ges_embedding(gestoken)

        epi = self.epi_enc(all_emb)

        x = torch.cat((tokens,positions,ges,epi),axis=2)

        x = self.do(x)
        #pdb.set_trace()
        
        x = self.tblocks(x)

        x = x.max(dim=1)[0] if self.max_pool else x.mean(dim=1) # pool over the time dimension

        feature = self.tofeature(x)

        if vis:
            return feature
        else:
            logits = self.toprobs(feature)

        # if we are given some desired targets also calculate the loss
        loss = None
        if targets is not None:
            loss = F.cross_entropy(logits.view(-1, logits.size(-1)), targets.view(-1))

        return logits, loss

class MuAtEpiF_EpiCompress(nn.Module):
    """
    Transformer for classifying sequences
    """
    def __init__(self, config):
        super().__init__()

        self.device = 'cpu'
        if torch.cuda.is_available():
            self.device = torch.cuda.current_device()

        self.num_tokens, self.max_pool = config.vocab_size, False

        self.token_embedding = nn.Embedding(config.vocab_size, config.n_embd,padding_idx=0)

        self.epiblocks = nn.Embedding(config.epi_size, config.epi_emb,padding_idx=0)
        self.epi_emb=config.epi_emb

        self.epi_enc = nn.Linear((config.epi_emb * 145), config.n_embd)
        
        tblocks = []
        for i in range(config.n_layer):
            tblocks.append(
                TransformerBlock(emb=int(config.n_embd), heads=config.n_head, seq_length=config.block_size, mask=False, dropout=config.attn_pdrop))

        self.tblocks = nn.Sequential(*tblocks)

        self.tofeature = nn.Sequential(nn.Linear(int(config.n_embd), 24),
                                        nn.ReLU())

        self.toprobs = nn.Linear(24, config.num_class)

        self.do = nn.Dropout(config.embd_pdrop)

    def forward(self, x,targets=None,vis=None):
        
        epitoken = x[0]

        batch = epitoken.shape[0]

        all_emb = torch.zeros(batch,epitoken.shape[1], self.epi_emb*epitoken.shape[2]).to(self.device)
        
        start = 0
        end = self.epi_emb
        for i in range(epitoken.shape[2]):
            oneinput = epitoken[:,:,i]
            oneres = self.epiblocks.to(self.device)(oneinput)
            all_emb[:,:,start:end] = oneres
            start = start + self.epi_emb
            end = end + self.epi_emb

        epi = self.epi_enc(all_emb)

        x = epi

        x = self.do(x)
        #pdb.set_trace()
        
        x = self.tblocks(x)

        x = x.max(dim=1)[0] if self.max_pool else x.mean(dim=1) # pool over the time dimension

        feature = self.tofeature(x)

        if vis:
            return feature
        else:
            logits = self.toprobs(feature)

        # if we are given some desired targets also calculate the loss
        loss = None
        if targets is not None:
            loss = F.cross_entropy(logits.view(-1, logits.size(-1)), targets.view(-1))

        return logits, loss


class MuAtMotifF(nn.Module):
    """
    Transformer for classifying sequences
    """
    def __init__(self, config):
        super().__init__()

        self.num_tokens, self.max_pool = config.vocab_size, False

        self.token_embedding = nn.Embedding(config.vocab_size, config.n_embd,padding_idx=0)

        tblocks = []
        for i in range(config.n_layer):
            tblocks.append(
                TransformerBlock(emb=config.n_embd, heads=config.n_head, seq_length=config.block_size, mask=False, dropout=config.attn_pdrop))

        self.tblocks = nn.Sequential(*tblocks)

        self.tofeature = nn.Sequential(nn.Linear(int(config.n_embd), 24),
                                        nn.ReLU())

        self.toprobs = nn.Linear(24, config.num_class)

        self.do = nn.Dropout(config.embd_pdrop)

    def forward(self, x,targets=None,vis=None):

        x = x[0]
    
        tokens = self.token_embedding(x)

        b, t, e = tokens.size()
        
        x = self.do(tokens)

        x = self.tblocks(x)

        x = x.max(dim=1)[0] if self.max_pool else x.mean(dim=1) # pool over the time dimension

        feature = self.tofeature(x)

        if vis:
            return feature
        else:
            logits = self.toprobs(feature)

        # if we are given some desired targets also calculate the loss
        loss = None
        if targets is not None:
            loss = F.cross_entropy(logits.view(-1, logits.size(-1)), targets.view(-1))

        return logits, loss

class MuAtMotifPosition(nn.Module):
    """
    Transformer for classifying sequences
    """
    def __init__(self, config):
        super().__init__()

        self.num_tokens, self.max_pool = config.vocab_size, False

        self.token_embedding = nn.Embedding(config.vocab_size, config.n_embd,padding_idx=0)
        self.position_embedding = nn.Embedding(config.position_size, config.n_embd,padding_idx=0)

        tblocks = []
        for i in range(config.n_layer):
            tblocks.append(
                TransformerBlock(emb=int(config.n_embd + config.n_embd), heads=config.n_head, seq_length=config.block_size, mask=False, dropout=config.attn_pdrop))

        self.tblocks = nn.Sequential(*tblocks)

        self.toprobs = nn.Linear(int(config.n_embd + config.n_embd), config.num_class)

        self.do = nn.Dropout(config.embd_pdrop)

    def forward(self, x,targets=None,vis=None):

        triplettoken = x[0]
        #pdb.set_trace()
        postoken = x[1]
        
    
        tokens = self.token_embedding(triplettoken)
        positions = self.position_embedding(postoken)

        x = torch.cat((tokens,positions),axis=2)

        x = self.do(x)
        
        x = self.tblocks(x)

        x = x.max(dim=1)[0] if self.max_pool else x.mean(dim=1) # pool over the time dimension

        logits = self.toprobs(x)

        # if we are given some desired targets also calculate the loss
        loss = None
        if targets is not None:
            loss = F.cross_entropy(logits.view(-1, logits.size(-1)), targets.view(-1))

        return logits, loss

class MuAtMotifPositionF(nn.Module):
    """
    Transformer for classifying sequences
    """
    def __init__(self, config):
        super().__init__()

        self.num_tokens, self.max_pool = config.vocab_size, False

        self.token_embedding = nn.Embedding(config.vocab_size, config.n_embd,padding_idx=0)
        self.position_embedding = nn.Embedding(config.position_size, config.n_embd,padding_idx=0)

        #pdb.set_trace()

        tblocks = []
        for i in range(config.n_layer):
            tblocks.append(
                TransformerBlock(emb=int(config.n_embd + config.n_embd), heads=config.n_head, seq_length=config.block_size, mask=False, dropout=config.attn_pdrop))

        self.tblocks = nn.Sequential(*tblocks)

        self.tofeature = nn.Sequential(nn.Linear(int(config.n_embd + config.n_embd), 24),
                                        nn.ReLU())

        #pdb.set_trace()
        self.toprobs = nn.Linear(24, config.num_class)

        self.do = nn.Dropout(config.embd_pdrop)

    def forward(self, x,targets=None,vis=None,visatt=None):

        triplettoken = x[0]
        postoken = x[1]
    
        tokens = self.token_embedding(triplettoken)
        positions = self.position_embedding(postoken)

        x = torch.cat((tokens,positions),axis=2)

        x = self.do(x)
        
        if visatt:
            dot1 = self.tblocks[0].attention(x,vis=True)

            manual_tblock = self.tblocks[0](x)

            after_tblock =self.tblocks(x)



            '''
            for block in self.tblocks:

                dot = block.attention(x,vis=True)

                pdb.set_trace()
                
                print(name)
            '''

            return dot1
        
        else:
            x = self.tblocks(x)

        x = x.max(dim=1)[0] if self.max_pool else x.mean(dim=1) # pool over the time dimension

        feature = self.tofeature(x)

        if vis:
            return feature
        else:
            logits = self.toprobs(feature)

        # if we are given some desired targets also calculate the loss
        loss = None
        if targets is not None:
            loss = F.cross_entropy(logits.view(-1, logits.size(-1)), targets.view(-1))

        return logits, loss


class MuAtMotifPositionGES(nn.Module):
    """
    Transformer for classifying sequences
    """
    def __init__(self, config):
        super().__init__()

        self.num_tokens, self.max_pool = config.vocab_size, False

        self.token_embedding = nn.Embedding(config.vocab_size, config.n_embd,padding_idx=0)
        self.position_embedding = nn.Embedding(config.position_size, config.n_embd,padding_idx=0)
        self.ges_embedding = nn.Embedding(config.ges_size+1, 4,padding_idx=0)

        tblocks = []
        for i in range(config.n_layer):
            tblocks.append(
                TransformerBlock(emb=int(config.n_embd + config.n_embd + 4), heads=config.n_head, seq_length=config.block_size, mask=False, dropout=config.attn_pdrop))

        self.tblocks = nn.Sequential(*tblocks)

        self.toprobs = nn.Linear(int(config.n_embd + config.n_embd + 4), config.num_class)

        self.do = nn.Dropout(config.embd_pdrop)

    def forward(self, x,targets=None,vis=None):

        triplettoken = x[0]
        postoken = x[1]
        gestoken = x[2]
    
        tokens = self.token_embedding(triplettoken)
        positions = self.position_embedding(postoken)
        ges = self.ges_embedding(gestoken)

        x = torch.cat((tokens,positions,ges),axis=2)

        x = self.do(x)
        #pdb.set_trace()
        
        x = self.tblocks(x)

        x = x.max(dim=1)[0] if self.max_pool else x.mean(dim=1) # pool over the time dimension

        logits = self.toprobs(x)

        # if we are given some desired targets also calculate the loss
        loss = None
        if targets is not None:
            loss = F.cross_entropy(logits.view(-1, logits.size(-1)), targets.view(-1))

        return logits, loss

class MuAtMotifPositionGESF(nn.Module):
    """
    Transformer for classifying sequences
    """
    def __init__(self, config):
        super().__init__()

        self.num_tokens, self.max_pool = config.vocab_size, False

        self.token_embedding = nn.Embedding(config.vocab_size, config.n_embd,padding_idx=0)
        self.position_embedding = nn.Embedding(config.position_size, config.n_embd,padding_idx=0)
        self.ges_embedding = nn.Embedding(config.ges_size+1, 4,padding_idx=0)

        tblocks = []
        for i in range(config.n_layer):
            tblocks.append(
                TransformerBlock(emb=int(config.n_embd + config.n_embd + 4), heads=config.n_head, seq_length=config.block_size, mask=False, dropout=config.attn_pdrop))

        self.tblocks = nn.Sequential(*tblocks)

        self.tofeature = nn.Sequential(nn.Linear(int(config.n_embd + config.n_embd + 4), 24),
                                        nn.ReLU())

        self.toprobs = nn.Linear(24, config.num_class)

        self.do = nn.Dropout(config.embd_pdrop)

    def forward(self, x,targets=None,vis=None):

        triplettoken = x[0]
        postoken = x[1]
        gestoken = x[2]
    
        tokens = self.token_embedding(triplettoken)
        positions = self.position_embedding(postoken)
        ges = self.ges_embedding(gestoken)

        x = torch.cat((tokens,positions,ges),axis=2)

        x = self.do(x)
        #pdb.set_trace()
        
        x = self.tblocks(x)

        x = x.max(dim=1)[0] if self.max_pool else x.mean(dim=1) # pool over the time dimension

        feature = self.tofeature(x)

        if vis:
            return feature
        else:
            logits = self.toprobs(feature)

        # if we are given some desired targets also calculate the loss
        loss = None
        if targets is not None:
            loss = F.cross_entropy(logits.view(-1, logits.size(-1)), targets.view(-1))

        return logits, loss

class MuAtMotifPositionGESRT(nn.Module):
    """
    Transformer for classifying sequences
    """
    def __init__(self, config):
        super().__init__()

        self.num_tokens, self.max_pool = config.vocab_size, False

        self.token_embedding = nn.Embedding(config.vocab_size, config.n_embd,padding_idx=0)
        self.position_embedding = nn.Embedding(config.position_size, config.n_embd,padding_idx=0)
        self.ges_embedding = nn.Embedding(config.ges_size+1, 4,padding_idx=0)

        tblocks = []
        for i in range(config.n_layer):
            tblocks.append(
                TransformerBlock(emb=int(config.n_embd + config.n_embd + 4 + 1), heads=config.n_head, seq_length=config.block_size, mask=False, dropout=config.attn_pdrop))

        self.tblocks = nn.Sequential(*tblocks)

        self.toprobs = nn.Linear(int(config.n_embd + config.n_embd + 4 + 1), config.num_class)

        self.do = nn.Dropout(config.embd_pdrop)

    def forward(self, x,targets=None,vis=None):

        triplettoken = x[0]
        postoken = x[1]
        gestoken = x[2]
        rt = x[3]
    
        tokens = self.token_embedding(triplettoken)
        positions = self.position_embedding(postoken)
        ges = self.ges_embedding(gestoken)

        x = torch.cat((tokens,positions,ges,rt),axis=2)

        x = self.do(x)
        
        x = self.tblocks(x)

        x = x.max(dim=1)[0] if self.max_pool else x.mean(dim=1) # pool over the time dimension

        logits = self.toprobs(x)

        # if we are given some desired targets also calculate the loss
        loss = None
        if targets is not None:
            loss = F.cross_entropy(logits.view(-1, logits.size(-1)), targets.view(-1))

        return logits, loss

class TransformerBlock(nn.Module):
    def __init__(self, emb, heads, mask, seq_length, ff_hidden_mult=4, dropout=0.0):
        super().__init__()

        self.attention = SelfAttention(emb, heads=heads, mask=mask)
        self.mask = mask

        self.norm1 = nn.LayerNorm(emb)
        self.norm2 = nn.LayerNorm(emb)

        self.ff = nn.Sequential(
            nn.Linear(emb, ff_hidden_mult * emb),
            nn.ReLU(),
            nn.Linear(ff_hidden_mult * emb, emb)
        )

        self.do = nn.Dropout(dropout)

    def forward(self, x):

        attended = self.attention(x)

        x = self.norm1(attended + x)

        x = self.do(x)

        fedforward = self.ff(x)

        x = self.norm2(fedforward + x)

        x = self.do(x)

        return x

class SelfAttention(nn.Module):
    def __init__(self, emb, heads=8, mask=False):
        """
        :param emb:
        :param heads:
        :param mask:
        """

        super().__init__()

        self.emb = emb
        self.heads = heads
        self.mask = mask

        self.tokeys = nn.Linear(emb, emb * heads, bias=False)
        self.toqueries = nn.Linear(emb, emb * heads, bias=False)
        self.tovalues = nn.Linear(emb, emb * heads, bias=False)

        self.unifyheads = nn.Linear(heads * emb, emb)
    
    def contains_nan(self,tensor):
        return bool((tensor != tensor).sum() > 0)

    def forward(self, x,vis=False):

        b, t, e = x.size()
        h = self.heads
        assert e == self.emb, f'Input embedding dim ({e}) should match layer embedding dim ({self.emb})'

        keys    = self.tokeys(x).view(b, t, h, e)
        queries = self.toqueries(x).view(b, t, h, e)
        values  = self.tovalues(x) .view(b, t, h, e)

        # compute scaled dot-product self-attention

        # - fold heads into the batch dimension
        keys = keys.transpose(1, 2).contiguous().view(b * h, t, e)
        queries = queries.transpose(1, 2).contiguous().view(b * h, t, e)
        values = values.transpose(1, 2).contiguous().view(b * h, t, e)

        # - get dot product of queries and keys, and scale
        dot = torch.bmm(queries, keys.transpose(1, 2))
        #pdb.set_trace()
        dot = dot / math.sqrt(e) # dot contains b*h  t-by-t matrices with raw self-attention logits

        if vis:
            return dot

        assert dot.size() == (b*h, t, t), f'Matrix has size {dot.size()}, expected {(b*h, t, t)}.'

        if self.mask: # mask out the lower half of the dot matrix,including the diagonal
            mask_(dot, maskval=float('-inf'), mask_diagonal=False)

        dot = F.softmax(dot, dim=2) # dot now has row-wise self-attention probabilities        

        assert not self.contains_nan(dot[:, 1:, :]) # only the forst row may contain nan

        if self.mask == 'first':
            dot = dot.clone()
            dot[:, :1, :] = 0.0
            # - The first row of the first attention matrix is entirely masked out, so the softmax operation results
            #   in a division by zero. We set this row to zero by hand to get rid of the NaNs

        # apply the self attention to the values
        out = torch.bmm(dot, values).view(b, h, t, e)

        # swap h, t back, unify heads
        out = out.transpose(1, 2).contiguous().view(b, t, h * e)


        return self.unifyheads(out)

if __name__ == '__main__':

    pdb.set_trace()


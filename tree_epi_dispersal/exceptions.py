from typing import Union


class InvalidDispersalSetup(Exception):

    def __init__(self, model: str, ell:Union[list, float, int], msg='The form of dispersal_type entered is incorrect!'):
        self.model = model
        self.ell = ell
        self.msg = msg

    def __str__(self):
        return self.msg + f'\n Parsed : {self.ell}, for model : {self.model}'



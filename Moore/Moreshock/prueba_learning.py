#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import boto3


client=boto3.client('s3')
path='s3://ifood-data-architect-test-source/restaurant.csv'
df=pd.read

res=pd.read_csv("/media/andres/DISPOSITIVO/mestrado/Andressss/lowpress_drpl/Mosesres/Moses434/resdrpl434.csv")